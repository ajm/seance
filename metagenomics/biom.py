import sys
import datetime

class BiomFile(object) :
    def __init__(self) :
        self.samples = []
        self.metadata = {}
        self.otus = []
        self.data = []

    def add_sample(self, sample_name, sample_metadata=None) :
        self.samples.append(sample_name)
        if sample_metadata :
            self.metadata[sample_name] = sample_metadata

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
        tmp["Date"] = md.get("date")
        tmp["Lemur"] = md.get("lemur")
        tmp["Location"] = md.get("location")
        tmp["Eggs"] = str(md.get("eggs"))
        return "{ " + ", ".join(map(lambda x: "%s : %s" % (repr(x[0]).replace("'","\""), repr(x[1]).replace("'","\"")), tmp.items())) + " }"

    def write_to(self, filename) :
        f = open(filename, 'w')
        self.write(f)
        f.close()

    def write(self, f=sys.stdout) :
        print >> f, '{'

        print >> f, '\t"id":null,'
        print >> f, '\t"format": "Biological Observation Matrix 0.9.1-dev",'
        print >> f, '\t"format_url": "http://biom-format.org/documentation/format_versions/biom-1.0.html",'
        print >> f, '\t"type": "OTU table",'
        print >> f, '\t"generated_by": "lemurs shitting in traps",'
        print >> f, '\t"date": "' + datetime.date.today().isoformat() + '",'

        print >> f, '\t"rows":['
        for i in self.otus :
            print >> f, "\t\t{ \"id\":\"%s\", \"metadata\":null }%s" % (i, "," if i != self.otus[-1] else "")
        print >> f, '\t],'
        
        print >> f, '\t"columns":['
        for i in self.samples :
            print >> f, "\t\t{ \"id\":\"%s\", \"metadata\": %s }%s" % (i, self.sample_metadata(i), "," if i != self.samples[-1] else "")
        print >> f, '\t],'

        print >> f, '\t"matrix_type": "sparse",'
        print >> f, '\t"matrix_element_type": "int",'

        print >> f, "\t\"shape\": [%d, %d]," % (len(self.otus), len(self.samples))
        
        print >> f, '\t"data":['
        for i in self.data :
            print >> f, "\t\t[%d, %d, %d]%s" % (i[0], i[1], i[2], "," if i != self.data[-1] else "")
        print >> f, '\t]'

        print >> f, '}'

