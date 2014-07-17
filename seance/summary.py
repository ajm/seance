from os.path import exists
from collections import defaultdict

class Summary(object) :
    def __init__(self, filename) :
        self.delim = ','
        self.data = defaultdict(dict)

        if exists(filename) :
            self.read(filename)

    def read(self, filename) :
        with open(filename) as f :
            header = f.readline().strip().split(self.delim)

            for line in f :
                line = line.strip()
                if line :
                    tmp = line.split(self.delim)
                    for ind,i in enumerate(tmp[1:]) :
                        self.data[tmp[0]][header[ind+1]] = i

    def update(self, other) :
        self.data.update(other)

    def __calc_totals(self) :
        totals = defaultdict(int)

        for i in self.data.keys() :
            if i == 'Totals' :
                continue

            for j in self.data[i] :
                tmp = self.data[i][j]

                if tmp == 'NA' :
                    continue

                totals[j] += int(tmp) 

        self.data['Totals'] = totals

    def __join(self, d) :
        return self.delim.join([ str(i) for i in d ])

    def write(self, filename) :
        self.__calc_totals()
        fields = sorted(self.data['Totals'].keys(), key=lambda x : x in ('Unique','Accepted'))

        def extract_in_order(d, fields) :
            return [ d.get(f, 'NA') for f in fields ]

        with open(filename, 'w') as f :
            print >> f, self.delim.join(['Filename'] + fields)

            for i in sorted(self.data.keys(), key=lambda x : (x, x in ('Totals'))) :
                print >> f, self.__join([i] + extract_in_order(self.data[i], fields))


if __name__ == '__main__' :
    s = Summary('summary.csv')
    s.write('test.csv')

