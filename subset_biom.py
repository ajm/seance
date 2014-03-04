import sys
import json

samples = ["CATI2",
        "ELE2",
        "VUL2",
        "SPI2",
        "NATI2",
        "EQU2",
        "X241",
        "X147",
        "X308",
        "X341",
        "X251",
        "X266",
        "X276",
        "X231",
        "X110",
        "X294",
        "X119",
        "X186",
        "N301",
        "N302",
        "N303",
        "M305",
        "M306",
        "X147+X308_1:1",
        "X241+X251_4:1",
        "X294+X231_10:1",
        "X251+X276+X119_1:1:1",
        "X231+X186+X266_10:2:1",
        "X147+X308_1:1_2",
        "X241+X251_4:1_2",
        "X294+X231_10:1_2",
        "X251+X276+X119_1:1:1_2",
        "X231+X186+X266_10:2:1_2"]

def include_sample(name) :
    return name.split()[0] in samples

def main() :
    f = open(sys.argv[1])
    biom = json.loads(f.read())
    f.close()

    #print json.dumps(biom)
    #return 0

    tmp = {}
    for key in biom :
        if key not in ('rows','columns','data','shape') :
            tmp[key] = biom[key]

    tmp['rows'] = []
    row_conversion = {}

    for index,row in enumerate(biom['rows']) :
        #if "Strongyloides" in row['id'] :
        if 1 :
            row_conversion[index] = len(tmp['rows'])
            tmp['rows'].append(row)

    tmp['columns'] = []
    column_conversion = {}

    for index,column in enumerate(biom['columns']) :
        metadata = column['metadata']
        #if (metadata['Location'] in ('Campsite', 'Talatakely')) and \
        #   (metadata['Lemur'] not in ('Dog', 'Rattus')) :
        if include_sample(column['id']) :
            column_conversion[index] = len(tmp['columns'])
            tmp['columns'].append(column)

    tmp['shape'] = [len(tmp['rows']), len(tmp['columns'])]

    tmp['data'] = []
    for rindex,cindex,value in biom['data'] :
        if (rindex in row_conversion) and (cindex in column_conversion) :
            tmp['data'].append([row_conversion[rindex], column_conversion[cindex], value])

    print json.dumps(tmp, sort_keys=True, indent=4)

    return 0

if __name__ == '__main__' :
    sys.exit(main())

