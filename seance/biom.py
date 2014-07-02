import sys
import datetime
import json

class BiomFile(object) :
    def __init__(self) :
        self.samples = []
        self.metadata = {}
        self.otus = []
        self.data = []

    def set_samples(self, sample_list) :
        self.samples = [ s.description() for s in sample_list ]
        self.metadata = dict(zip(self.samples, [ s.metadata for s in sample_list ]))

    def add_sample(self, sample_name, sample_metadata=None) :
        self.samples.append(sample_name)
        if sample_metadata :
            self.metadata[sample_name] = sample_metadata

    def set_otus(self, otu_list) :
        self.otus = otu_list

    def add_otu(self, otu_name) :
        self.otus.append(otu_name)

    def add_quantity(self, otuid, sampleid, value) :
        if value > 0 :
            self.data.append((otuid, sampleid, value))

    def sample_metadata(self, sample_name) :
        if sample_name not in self.metadata :
            return "null"
        tmp = {}
        md = self.metadata[sample_name]
        tmp["Date"] = str(md["date"])
        tmp["Lemur"] = md["lemur"]
        tmp["Location"] = md["location"]
        tmp["Eggs"] = str(md["eggs"])
        return tmp

    def sample_metadata_str(self, sample_name) :
        tmp = self.sample_metadata(sample_name)
        return "{ " + ", ".join(map(lambda x: "%s : %s" % (repr(x[0]).replace("'","\""), repr(x[1]).replace("'","\"")), tmp.items())) + " }"

    def write_to(self, filename) :
        f = open(filename, 'w')
        self.write(f)
        f.close()

    def write(self, f=sys.stdout) :
        tmp = {}

        tmp["id"] = "null"
        tmp["format"] = "Biological Observation Matrix 0.9.1-dev"
        tmp["format_url"] = "http://biom-format.org/documentation/format_versions/biom-1.0.html"
        tmp["type"] = "OTU table"
        tmp["generated_by"] = "seance"
        tmp["date"] = datetime.date.today().isoformat()
        tmp["rows"] = []
        tmp["columns"] = []
        tmp["matrix_type"] = "sparse"
        tmp["matrix_element"] = "int"
        tmp["shape"] = [len(self.otus), len(self.samples)]
        tmp["data"] = []

        for id,label in self.otus :
            tmp["rows"].append({ "id" : id, "metadata" : { "label" : label } })

        for i in self.samples :
            tmp["columns"].append({ "id" : i, "metadata" : self.sample_metadata(i) })

        for i in self.data :
            tmp["data"].append(i)


        print >> f, json.dumps(tmp, sort_keys=True, indent=2, separators=(",",": "))

    def change_otu_names(self, filename, names) :
        tmp = json.load(open(filename))

        # iterate through and alter any in 'names'
        for r in tmp['rows'] :
            if r['id'] in names :
                r['metadata']['label'] = names[r['id']]

        # rewrite the file
        f = open(filename, 'w')
        print >> f, json.dumps(tmp, sort_keys=True, indent=2, separators=(",",": "))
        f.close()

    def get_label_mapping(self, filename) :
        data = json.load(open(filename))
        tmp = {}

        for r in data['rows'] :
            tmp[r['id']] = r['metadata']['label']

        return tmp

