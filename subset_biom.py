import sys
import json

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
        if "Strongyloides" in row['id'] :
            row_conversion[index] = len(tmp['rows'])
            tmp['rows'].append(row)

    tmp['columns'] = []
    column_conversion = {}

    for index,column in enumerate(biom['columns']) :
        metadata = column['metadata']
        if (metadata['Location'] in ('Campsite', 'Talatakely')) and \
           (metadata['Lemur'] not in ('Dog', 'Rattus')) :
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

