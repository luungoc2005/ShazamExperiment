using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ShazamExperiment.Models
{
    public class SongData
    {
        [JsonProperty("trackName")]
        public string TrackName { get; set; }

        [JsonProperty("fileName")]
        public string FileName { get; set; }

        [JsonProperty("songFingerprint")]
        public IEnumerable<long> SongFingerprint { get; set; }

        public void SaveToFile(string fileName)
        {
            if (File.Exists(fileName)) File.Delete(fileName);
            var jsonString = JsonConvert.SerializeObject(this);
            using (var writer = new StreamWriter(File.Open(fileName, FileMode.CreateNew)))
            {
                writer.Write(jsonString);
                writer.Flush();
            }
        }
    }
}
