using System;
using System.Collections.Generic;
using System.Windows.Forms;
using NAudio.Wave;
using System.IO;
using NAudio.Utils;
using System.Numerics;
using ShazamExperiment.Models;
using Newtonsoft.Json;
using System.Linq;
using System.Drawing;

namespace ShazamExperiment
{
    public partial class Form1 : Form
    {
        private struct SongTimeData
        {
            public int Position;
            public SongData SongData;
        }
        private List<SongData> songList = new List<SongData>();
        private Dictionary<long, IEnumerable<SongTimeData>> songDb = new Dictionary<long, IEnumerable<SongTimeData>>();

        DirectoryInfo _currentDirectory;
        string CurrentDirectory
        {
            get
            {
                return (_currentDirectory ?? new DirectoryInfo(Directory.GetCurrentDirectory())).FullName;
            }
            set
            {
                if (Directory.Exists(value))
                {
                    _currentDirectory = new DirectoryInfo(value);
                }

                //update everything
                label1.Text = CurrentDirectory;
                var itemList = new List<ListViewItem>();
                foreach (var file in (new DirectoryInfo(CurrentDirectory)).GetFiles())
                {
                    if (!IsSongExist(file.Name) && Path.GetExtension(file.Name).ToLower() == ".mp3")
                    {
                        var listItem = new ListViewItem(file.Name);
                        listItem.SubItems.Add("");
                        itemList.Add(listItem);
                    }
                }
                listView1.Items.AddRange(itemList.ToArray());
                UpdateList();
            }
        }

        string _dataDirectory = "data";
        DirectoryInfo DataDirectory
        {
            get
            {
                var dataPath = Path.Combine(Directory.GetCurrentDirectory(), _dataDirectory);
                if (!Directory.Exists(dataPath))
                {
                    return Directory.CreateDirectory(dataPath);
                }
                else return new DirectoryInfo(dataPath);
            }
        }

        private string GetDataFilePath(string fileName)
        {
            return Path.Combine(DataDirectory.FullName, fileName);
        }

        //private static int LOWER_LIMIT = 0;
        private static int UPPER_LIMIT = 1024;
        private static short DAMPENING_FACTOR = 2;
        private static int[] RANGE = new int[] { 10, 20, 40, 80, 160, 511, UPPER_LIMIT + 1 };

        private WaveIn RecordStream { get; set; }
        private MemoryStream RecordMemoryStream { get; set; }

        private System.Diagnostics.Stopwatch Stopwatch { get; set; } = new System.Diagnostics.Stopwatch();

        private Complex[][] SpectrumData { get; set; }

        public Form1()
        {
            InitializeComponent();
        }

        private bool IsSongExist(string fileName)
        {
            return songList.Where(x => (x.TrackName == fileName)).Any();
        }

        private void UpdateList()
        {
            foreach(ListViewItem item in listView1.Items)
            {
                if (IsSongExist(item.SubItems[0].Text))
                {
                    item.SubItems[1].Text = "Yes";
                }
                else
                {
                    item.SubItems[1].Text = "No";
                }
            }
        }

        //private static double ComplexAbsolute(Complex input)
        //{
        //    try
        //    {
        //        return Math.Sqrt(input.X * input.X + input.Y * input.Y);
        //        //if (input.Y == 0)
        //        //    return 0;
        //        //else if (input.X > input.Y)
        //        //    return input.X * Math.Sqrt(1 + (input.Y * input.Y) / (input.X * input.X));
        //        //else
        //        //    return input.Y * Math.Sqrt(1 + (input.X * input.X) / (input.Y * input.Y));
        //    }
        //    catch
        //    {
        //        return double.NaN;
        //    }
        //}

        private SongData AddDataFile(string fileName)
        {
            if (!File.Exists(fileName)) return null;
            var stringData = string.Empty;
            using (var reader = new StreamReader(File.Open(fileName, FileMode.Open)))
            {
                stringData = reader.ReadToEnd();
            }
            try
            {
                var returnObject = JsonConvert.DeserializeObject<SongData>(stringData);
                if (returnObject != null)
                {
                    AddData(returnObject);
                }
                songList.Add(returnObject);
                return returnObject;
            }
            catch
            {
                return null;
            }
        }

        private void AddData(SongData data)
        {
            if (data.SongFingerprint != null && data.SongFingerprint.Any())
            {
                int count = 0;
                foreach (var fingerprint in data.SongFingerprint)
                {
                    IEnumerable<SongTimeData> dataList;
                    if (!songDb.TryGetValue(fingerprint, out dataList))
                    {
                        dataList = new List<SongTimeData>();
                        songDb.Add(fingerprint, dataList);
                    }
                    ((List<SongTimeData>)dataList).Add(new SongTimeData()
                    {
                        SongData = data,
                        Position = count
                    });
                    count++;
                }
            }
        }

        private void button3_Click(object sender, EventArgs e)
        {
            var fileList = new DirectoryInfo(CurrentDirectory).GetFiles();
            progressBar1.Maximum = fileList.Length;
            progressBar1.Value = 0;
            foreach (var file in fileList)
            {
                if (Path.GetExtension(file.FullName).ToLower() == ".mp3")
                {
                    if (!IsSongExist(file.Name))
                    {
                        label2.Text = $"Analyzing {file.Name} - {progressBar1.Value}/{progressBar1.Maximum}";
                        Application.DoEvents();
                        var newSong = new SongData()
                        {
                            TrackName = file.Name,
                            FileName = file.FullName,
                            SongFingerprint = GetFileFingerprint(file.FullName)
                        };
                        newSong.SaveToFile(GetDataFilePath(file.Name + ".txt"));
                        AddData(newSong);
                    }
                }
                progressBar1.Value += 1;
                //UpdateList();
                Application.DoEvents();
            }
        }

        private long GetHash(Array data)
        {
            if (data.Length < 2) return 0;
            long retVal = 0;
            bool doReturn = false;
            for (int i = data.Length - 2; i > -1; i--)
            {
                var point = Convert.ToInt64(data.GetValue(i));
                retVal *= 1000;
                retVal += point;
                if (retVal == 0) retVal = 1;
                if (point != 0) doReturn = true;
            }
            return doReturn ? retVal : 0;
        }

        private List<long> GetFingerprint(IWaveProvider samples)
        {
            //write to file for testing
            //using (var file = File.Open(Guid.NewGuid() + ".wav", FileMode.CreateNew))
            //    {
            //        using (var writer = new WaveFileWriter(file, samples.WaveFormat))
            //        {
            //            var fileBuffer = new float[4096];
            //            var bufferLength = 0;
            //            while ((bufferLength = samples.Read(fileBuffer, 0, fileBuffer.Length)) > 0)
            //            {
            //                writer.WriteSamples(fileBuffer, 0, bufferLength);
            //            }
            //        }
            //    }
            var retList = new List<long>();
            var getIndex = new Func<int, int>(x =>
            {
                int i = 0;
                while (RANGE[i] < x) i++;
                return i;
            });
            int windowLength = 4096, length = 0;
            var buffer = new byte[windowLength * 4];
            var wrapper = new WaveBuffer(buffer);
            var floatBuffer = new float[windowLength];
            var results = new List<Complex[]>();
            int m = (int)Math.Log(windowLength, 2.0);
            while ((length = samples.Read(buffer, 0, buffer.Length)) > 0)
            {
                //Buffer.BlockCopy(buffer, 0, floatBuffer, 0, length);
                var complex = new Complex[windowLength];
                for (int i = 0; i < length /4; i++)
                {
                    //complex[i] = new Complex(floatBuffer[i], 0);
                    complex[i] = new Complex(wrapper.FloatBuffer[i], 0);
                }
                //FastFourierTransform.FFT(true, m, complex);
                Accord.Math.FourierTransform.FFT(complex, Accord.Math.FourierTransform.Direction.Forward);
                results.Add(complex);
            }
            foreach (var result in results)
            {
                if (result != null && result.Any())
                {
                    //System.Diagnostics.Debug.WriteLine(result.Length);
                    var highscores = new double[RANGE.Length];
                    var recordPoints = new double[RANGE.Length];
                    for (int freq = 0; freq < Math.Min(UPPER_LIMIT + 1, result.Length); freq++)
                    {
                        double mag = Math.Log(result[freq].Magnitude + 1);
                        int index = getIndex(freq);
                        if (mag > highscores[index])
                        {
                            highscores[index] = mag;
                            recordPoints[index] = Math.Max(0, freq - (freq % DAMPENING_FACTOR));
                        }
                    }
                    var hash = GetHash(recordPoints);
                    if (hash != 0) retList.Add(GetHash(recordPoints));
                }
            }
            return retList;
        }

        private List<long> GetFingerprint(Stream inputStream) //Assumes mono channel, 32 bit
        {
            inputStream.Seek(0, SeekOrigin.Begin);
            var retList = new List<long>();
            var getIndex = new Func<int, int>(x =>
            {
                int i = 0;
                while (RANGE[i] < x) i++;
                return i;
            });
            var isPowerOfTwo = new Func<int, bool>(x => ((x & (x - 1)) == 0));
            int length = 0;
            int windowLength = 2048;
            int windowBytesLength = windowLength * 4;
            var buffer = new byte[windowBytesLength];
            var floatBuffer = new float[windowLength];
            var results = new Complex[inputStream.Length / windowBytesLength + 1][];
            int m = (int)Math.Log(windowLength, 2.0);
            int times = 0;
            while ((length = inputStream.Read(buffer, 0, buffer.Length)) > 0)
            {
                var wrapper = new WaveBuffer(buffer);
                var complex = new Complex[windowLength];
                //Buffer.BlockCopy(buffer, 0, floatBuffer, 0, length);
                for (int i = 0; i < length / 4; i++)
                {
                    complex[i] = new Complex(wrapper.FloatBuffer[i] * (float)NAudio.Dsp.FastFourierTransform.HammingWindow(i, windowLength), 0);
                    //complex[i].Y = 0;
                    //complex[i] = new Complex(floatBuffer[i], 0);
                    //System.Diagnostics.Debug.Print(wrapper.FloatBuffer[i].ToString());
                }
                Accord.Math.Transforms.FourierTransform2.FFT(complex, Accord.Math.FourierTransform.Direction.Forward);
                //for (int i = 0; i < length; i++)
                //{
                //    switch (sampleSize)
                //    {
                //        case 16:
                //            complex[i].X = wrapper.ShortBuffer[i] * (float)FastFourierTransform.HammingWindow(i, windowLength);
                //            break;
                //        case 32:
                //            complex[i].X = wrapper.FloatBuffer[i] * (float)FastFourierTransform.HammingWindow(i, windowLength);
                //            break;
                //        case 64:
                //            complex[i].X = wrapper.IntBuffer[i] * (float)FastFourierTransform.HammingWindow(i, windowLength);
                //            break;
                //        default:
                //            complex[i].X = wrapper.ByteBuffer[i] * (float)FastFourierTransform.HammingWindow(i, windowLength);
                //            break;
                //    }
                //    complex[i].Y = 0;
                //}
                //FastFourierTransform.FFT(true, m, complex);
                results[times] = complex;
                times++;
            }
            //SpectrumData = results; // draw spectrum data
            //DrawSpectrum(Guid.NewGuid() + ".png");
            foreach (var result in results)
            {
                if (result != null && result.Any())
                {
                    //System.Diagnostics.Debug.WriteLine(result.Length);
                    var highscores = new double[RANGE.Length];
                    var recordPoints = new double[RANGE.Length];
                    for (int freq = 0; freq < UPPER_LIMIT + 1; freq++)
                    {
                        double mag = Math.Log(result[freq].Magnitude + 1);
                        int index = getIndex(freq);
                        if (mag > highscores[index])
                        {
                            highscores[index] = mag;
                            recordPoints[index] = Math.Max(0, freq - (freq % DAMPENING_FACTOR));
                        }
                    }
                    var hash = GetHash(recordPoints);
                    if (hash != 0) retList.Add(GetHash(recordPoints));
                }
            }
            GC.Collect();
            return retList;
        }

        private void CombineDualChannel(byte[] buffer, Stream stream, int length = -1, int sampleSize = 16)
        {
            var tempBuffer = new byte[0];
            if (length == -1) length = buffer.Length;

            var _sourceBuffer = new WaveBuffer(buffer);

            tempBuffer = BufferHelpers.Ensure(tempBuffer, length / 2);
            var _destBuffer = new WaveBuffer(tempBuffer);
            int pos = 0;
            int bytesSampleSize = sampleSize / 8;

            switch (sampleSize)
            {
                case 16:
                    for (int i = 0; i < length / bytesSampleSize; i += 2)
                    {
                        var outSample = _sourceBuffer.ShortBuffer[i] * 0.5f +
                                        _sourceBuffer.ShortBuffer[i + 1] * 0.5f;
                        outSample = Math.Max(Int16.MinValue, Math.Min(outSample, Int16.MaxValue));
                        _destBuffer.ShortBuffer[pos] = Convert.ToInt16(outSample);
                        pos++;
                    }
                    break;
                case 32:
                    for (int i = 0; i < length / bytesSampleSize; i += 2)
                    {
                        var left = _sourceBuffer.FloatBuffer[i];
                        var right = _sourceBuffer.FloatBuffer[i];
                        var outSample = _sourceBuffer.FloatBuffer[i] * 0.5f +
                                    _sourceBuffer.FloatBuffer[i + 1] * 0.5f;
                        _destBuffer.FloatBuffer[pos] = outSample;
                        pos++;
                    }
                    //outSample = Math.Max(float.MinValue, Math.Min(outSample, float.MaxValue));
                    break;
                case 64:
                    for (int i = 0; i < length / bytesSampleSize; i += 2)
                    {
                        var outSample = (_sourceBuffer.IntBuffer[i] * 0.5 +
                                        _sourceBuffer.IntBuffer[i + 1] * 0.5);
                        outSample = Math.Max(Int32.MinValue, Math.Min(outSample, Int32.MaxValue));
                        _destBuffer.IntBuffer[pos] = Convert.ToInt32(outSample);
                        pos++;
                    }
                    break;
                default:
                    for (int i = 0; i < length / bytesSampleSize; i += 2)
                    {
                        var outSample = _sourceBuffer.ByteBuffer[i] * 0.5f +
                                    _sourceBuffer.ByteBuffer[i + 1] * 0.5f;
                        outSample = Math.Max(byte.MinValue, Math.Min(outSample, byte.MaxValue));
                        _destBuffer.ByteBuffer[pos] = Convert.ToByte(outSample);
                        pos++;
                    }
                    break;
            }

            stream.Write(tempBuffer, 0, length / 2);
        }

        private List<long> GetFileFingerprint(string fileName)
        {
            if (!File.Exists(fileName)) return null;

            using (var memStream = new MemoryStream())
            {
                var buffer = new byte[65536];
                //var shortBuffer = new short[32768];
                var floatBuffer = new byte[65536];
                var tempBuffer = new byte[0];
                int length = 0;

                var floatWrapper = new WaveBuffer(floatBuffer);
                
                using (var reader = new Mp3FileReader(fileName))
                {
                    var bitsPerSample = reader.WaveFormat.BitsPerSample;

                    System.Diagnostics.Debug.WriteLine($"Analyzing {fileName} - {reader.WaveFormat.ToString()}");
                    //var desiredFormat = new WaveFormat(44100, 32, 1);
                    //var pcmStream = (reader.WaveFormat.Encoding != WaveFormatEncoding.Pcm) ?
                    //    WaveFormatConversionStream.CreatePcmStream(reader) : reader;
                    //var monoStream = new WaveFormatConversionStream(desiredFormat, pcmStream).ToSampleProvider(); 
                    //var provider = new StereoToMonoProvider16(reader);                   
                    //return GetFingerprint(provider);
                    //while ((length = monoStream.Read(floatBuffer, 0, floatBuffer.Length)) > 0)
                    //{
                    //    Buffer.BlockCopy(floatBuffer, 0, buffer, 0, length * 4);
                    //    memStream.Write(buffer, 0, length * 4);
                    //}
                    //MessageBox.Show(reader.WaveFormat.Encoding.ToString() + "-" + reader.WaveFormat.BitsPerSample);
                    //while ((length = reader.Read(buffer, 0, buffer.Length)) > 0)
                    //{
                    //    if (reader.WaveFormat.Channels == 2)
                    //    {
                    //        CombineDualChannel(buffer, memStream, length, bitsPerSample);
                    //    }
                    //    else
                    //    {
                    //        memStream.Write(buffer, 0, length);
                    //    }
                    //}
                    while ((length = reader.Read(buffer, 0, buffer.Length)) > 0)
                    {
                        //Buffer.BlockCopy(buffer, 0, shortBuffer, 0, length);
                        var wrapper = new WaveBuffer(buffer);
                        var count = 0;
                        for (int i = 0; i < length / 2; i += 2) //short buffer length will equal half of byte buffer
                        {
                            double outSample = wrapper.ShortBuffer[i] * 0.5 + wrapper.ShortBuffer[i + 1] * 0.5;
                            outSample = Math.Max(outSample, float.MinValue);
                            outSample = Math.Min(outSample, float.MaxValue);
                            floatWrapper.FloatBuffer[count++] = (float)outSample;
                            //memStream.Write(BitConverter.GetBytes((float)outSample), 0, sizeof(float));
                            //writer.Write((float)outSample);
                            //System.Diagnostics.Debug.WriteLine(outSample);
                        }
                        memStream.Write(floatWrapper.ByteBuffer, 0, length);
                    }
                    return GetFingerprint(memStream);
                }
            }
        }

        private void button4_Click(object sender, EventArgs e)
        {
            folderBrowserDialog1.SelectedPath = CurrentDirectory;
            folderBrowserDialog1.ShowDialog();
            CurrentDirectory = folderBrowserDialog1.SelectedPath;
        }

        private void Form1_Load(object sender, EventArgs e)
        {
            //load data
            listView1.BeginUpdate();
            foreach (var file in DataDirectory.GetFiles())
            {
                var fileAdded = AddDataFile(file.FullName);
                var listItem = new ListViewItem(fileAdded.TrackName);
                listItem.SubItems.Add("");
                listView1.Items.Add(listItem);
            }
            listView1.EndUpdate();
            UpdateList();

            int count = WaveIn.DeviceCount;
            for (int device = 0; device < count; device++)
            {
                WaveInCapabilities deviceInfo = WaveIn.GetCapabilities(device);
                comboBox1.Items.Add(deviceInfo.ProductName);
            }

            textBox1.Text = "test\r\ntest";
        }

        private void button1_Click(object sender, EventArgs e)
        {
            if (comboBox1.SelectedIndex > -1)
            {
                RecordStream = new WaveIn();
                RecordStream.DeviceNumber = comboBox1.SelectedIndex;
                WaveInCapabilities deviceInfo = WaveIn.GetCapabilities(RecordStream.DeviceNumber);
                RecordStream.WaveFormat = new WaveFormat(44100, 16, 2);

                if (RecordMemoryStream != null)
                {
                    RecordMemoryStream.Close();
                }
                RecordMemoryStream = new MemoryStream();

                RecordStream.DataAvailable += RecordStream_DataAvailable;
                RecordStream.StartRecording();

                textBox1.Text = "Started recording with " + comboBox1.SelectedItem.ToString();

                Stopwatch.Restart();
            }
        }

        private void RecordStream_DataAvailable(object sender, WaveInEventArgs e)
        {
            if (RecordMemoryStream != null)
            {
                RecordMemoryStream.Write(e.Buffer, 0, e.BytesRecorded);
                //CombineDualChannel(e.Buffer, RecordMemoryStream, e.Buffer.Length);
            }
        }

        private void button2_Click(object sender, EventArgs e)
        {
            if (RecordStream != null)
            {
                RecordStream.StopRecording();
                Stopwatch.Stop();

                if (RecordMemoryStream.Length > 0)
                {
                    RecordMemoryStream.Seek(0, SeekOrigin.Begin);
                    var stream = new RawSourceWaveStream(RecordMemoryStream, RecordStream.WaveFormat);

                    // find the song
                    textBox1.Text = $"Fingerprinting {RecordMemoryStream.Length} bytes, format: {RecordStream.WaveFormat.ToString()}... \r\n";
                    Application.DoEvents();
                    //var returnList = GetFingerprint(RecordMemoryStream, RecordStream.WaveFormat.BitsPerSample);

                    List<long> returnList;
                    using (var tempStream = new MemoryStream())
                    {
                        var buffer = new byte[65536];
                        var shortBuffer = new short[32768];
                        var length = 0;
                        while ((length = stream.Read(buffer, 0, buffer.Length)) > 0)
                        {
                            Buffer.BlockCopy(buffer, 0, shortBuffer, 0, length);
                            for (int i = 0; i < length / 2; i += 2)
                            {
                                double outSample = shortBuffer[i] * 0.5 + shortBuffer[i + 1] * 0.5;
                                outSample = Math.Max(outSample, float.MinValue);
                                outSample = Math.Min(outSample, float.MaxValue);
                                tempStream.Write(BitConverter.GetBytes((float)outSample), 0, sizeof(float));
                                //writer.Write((float)outSample);
                                //System.Diagnostics.Debug.WriteLine(outSample);
                            }
                        }
                        returnList = GetFingerprint(tempStream);
                    }

                    textBox1.Text += "Searching... \r\n";
                    Application.DoEvents();
                    var foundList = new List<KeyValuePair<SongData, int>>();
                    var noMatches = new List<KeyValuePair<SongData, int>>();
                    var tempList = new List<SongTimeData>();
                    foundList.Clear();
                    tempList.Clear();

                    int fingerprintPos = 0;
                    foreach (var fingerprint in returnList)
                    {
                        IEnumerable<SongTimeData> matchList;
                        if (songDb.TryGetValue(fingerprint, out matchList))
                        {
                            tempList.AddRange(matchList.Select(x =>
                            {
                                x.Position -= fingerprintPos;
                                return x;
                            }));
                        }
                        fingerprintPos++;
                    }

                    textBox1.Text += $"Found {tempList.Count} matches...\r\n";
                    Application.DoEvents();
                    foreach (var item in tempList)
                    {
                        if (!foundList.Where(x => x.Key.FileName ==item.SongData.FileName).Any())
                        {
                            var foundItem = tempList.Where(x => (x.SongData.FileName == item.SongData.FileName)).ToArray();
                            if (foundItem != null && foundItem.Any())
                            {
                                var minIndex = foundItem.Min(x => x.Position);
                                int score = 0;
                                int matches = foundItem.Length;
                                for (int i = 0; i < foundItem.Length; i++)
                                {
                                    score += Math.Abs(foundItem[i].Position - minIndex);
                                }
                                foundList.Add(new KeyValuePair<SongData, int>(item.SongData, score / matches));
                                noMatches.Add(new KeyValuePair<SongData, int>(item.SongData, matches));
                            }
                        }
                    }

                    SongData bestMatch = null;
                    int minScore = int.MaxValue;
                    int maxMatches = 0;
                    int count = 0;

                    foreach (var item in noMatches.OrderByDescending(x => x.Value))
                    {
                        var score = foundList.Where(x => (item.Key.FileName == x.Key.FileName)).FirstOrDefault();
                        if ((item.Value >= 2)) {
                            
                            if (count < 20)
                            {
                                textBox1.Text += $"{count}: Match {item.Key.TrackName} {item.Value} times - Score: {score.Value}\r\n";
                            }
                            count++;

                            if (item.Value > maxMatches || ((item.Value == maxMatches) && score.Value < minScore))
                            {
                                bestMatch = item.Key;
                                minScore = score.Value;
                                maxMatches = item.Value;
                            }
                        }                                
                    }
                    if (bestMatch != null) textBox1.Text += $"\r\nBest match: {bestMatch.TrackName}";
                }
            }
        }

        private void timer1_Tick(object sender, EventArgs e)
        {
            label3.Text = Stopwatch.Elapsed.ToString();
        }

        private void DrawSpectrum(Complex[][] spectrumData, string fileName)
        {
            if (spectrumData == null || !spectrumData.Any()) return;
            var lengthY = 2;
            var lengthX = 2;
            var newBitmap = new Bitmap(
                lengthX * spectrumData[0].Length / 2, 
                lengthY * spectrumData.Length, 
                System.Drawing.Imaging.PixelFormat.Format24bppRgb);
            using (var g = Graphics.FromImage(newBitmap))
            {
                g.CompositingQuality = System.Drawing.Drawing2D.CompositingQuality.HighQuality;
                g.CompositingMode = System.Drawing.Drawing2D.CompositingMode.SourceCopy;
                g.InterpolationMode = System.Drawing.Drawing2D.InterpolationMode.HighQualityBicubic;
                
                //g.Clear(Color.Green);
                for (int y = 0; y < SpectrumData.Length; y++)
                {
                    if (SpectrumData[y] != null && SpectrumData[y].Any())
                    {
                        for (int x = 0; x < SpectrumData[y].Length / 2; x++)
                        {
                            var abs = SpectrumData[y][x].Magnitude;
                            var mag = Math.Log(SpectrumData[y][x].Magnitude + 1);
                            g.FillRectangle(new SolidBrush(
                                    Color.FromArgb(0, Math.Min((int)mag * 10, 255), Math.Min((int)mag * 20, 255))),
                                x * lengthX, y * lengthY, 1, lengthY);
                        }
                    }
                }
            }
            newBitmap.Save(fileName, System.Drawing.Imaging.ImageFormat.Png);
            newBitmap.Dispose();
        }

        private void pictureBox1_Paint(object sender, PaintEventArgs e)
        {
            if (SpectrumData == null || !SpectrumData.Any()) return;

        }
    }
}
