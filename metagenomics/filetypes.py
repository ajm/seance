import sys
import os
from metagenomics.datatypes import Sequence, SampleMetaData, IUPAC

class DataFileError(Exception):
    pass

class DataFile(object) :
    def __init__(self, fname, extension) :
        self.pathname = os.path.dirname(fname)
        self.filename = os.path.basename(fname)
        self.extension = extension

        if not self.__exists() :
            raise DataFileError("'%s' does not exist" % self.get_filename())

    def __exists(self) :
        return os.path.isfile(self.get_filename())

    def get_filename(self) :
        return self.pathname + os.sep + self.filename

    def get_basename(self) :
        return self.filename

    def get_extension(self) :
        return self.extension

class SffFile(DataFile) :
    def __init__(self, fname) :
        DataFile.__init__(self, fname, ".sff")

class State(object) :
    def __init__(self, states) :
        self.__counter = 0
        self.__num_states = states

    def inc(self) :
        self.__counter += 1
        self.__counter %= self.__num_states

    def get(self) :
        return self.__counter

class FastqFile(DataFile) :
    def __init__(self, fname) :
        super(FastqFile, self).__init__(fname, ".fastq")
        self._filehandle = None
        self._sequences = {}

        self._validators = {
                    0 : self.__validate_seqid,
                    1 : self.__validate_sequence,
                    2 : self.__validate_qualid,
                    3 : self.__validate_qualities
                }

        self._state = State(4)
        self._linenum = 0
        self._current = {  0 : None,
                           1 : None,
                           2 : None,
                           3 : None }

    def __validate_seqid(self, s) :
        if not s.startswith('@') :
            raise ParseError("expected line %d to start with an @" % self._linenum)

    def __validate_sequence(self, s) :
        uniq = set(s)

        for i in uniq :
            if i not in IUPAC.codes :
                raise ParseError("line %d contained an invalid UIPAC code (%s)" % (self._linenum, i))

    def __validate_qualid(self, s) :
        if not s.startswith('+') :
            raise ParseError("expected line %d to start with an +" % self._linenum)

    def __validate_qualities(self, s) :
        uniq = set(s)

        for i in uniq :
            try :
                Sequence.convert_quality(i)
        
            except ValueError, ve :
                raise ParseError("line %d contained an invalid quality value (%s)" % (self._linenum, i))

    def __validate(self, s, state) :
        self._validators[state](s)

    def __iter__(self) :
        return self

    def next(self) :
        return self.read()

    def open(self) :
        self._filehandle = open(self.get_filename())

    def close(self) :
        self._filehandle.close()

    def read(self) :
        #try :
        #    f = open(self.__filename)

        #except IOError, ioe :
        #    raise ParseError(str(ioe))

        #self.__linenum = 1
        #st = State(4)
        
        #current = { 0 : None,
        #            1 : None,
        #            2 : None,
        #            3 : None }
        
        #for line in f :
        for line in self._filehandle :
            line = line.strip()
            
            self.__validate(line, self._state.get())
            self._current[self._state.get()] = line

            self._state.inc()
            self._linenum += 1

            #if self._state.get() == 3 :
            if self._state.get() == 0 :
                #self.__sequences[current[0][1:]] = Sequence(current[1], current[3])
                return Sequence(self._current[1], self._current[3])

            #self._state.inc()
            #self._linenum += 1

        #f.close()

        raise StopIteration

class MetaDataReader(object) :
    def __init__(self, metadata_fname) :
        self.metadata_fname = metadata_fname
        self.metadata = {}

    def get(self, key) :
        return self.metadata.get(key, None)

    def process(self) :
        num_required_fields = 7

        line_num = 0

        f = open(self.metadata_fname) 

        for line in f :
            line_num += 1
            line = line.strip()
            if line == '' :
                continue

            data = line.split()
            if len(data) != num_required_fields :
                print >> sys.stderr, \
                         "Warning: line %d of metadata file '%s' contains %d fields, expected %d" % \
                         (line_num, self.metadata_fname, len(data), num_required_fields)
                continue

            # Tg_2 53 Tg_2-25062012-1 Solofo 10/11/2010 Campsite 9
            filename = data[2] + ".sff"
            lemurname = data[3]
            date = data[4]
            location = data[5]
            try :
                num_eggs = int(data[6])

            except ValueError, ve :
                print >> sys.stderr, \
                         "Warning: line %d of metadata file '%s' egg field is not a number (read '%d')" % \
                         (line_num, self.metadata_file, data[6])
                continue
            
            smd = SampleMetaData()
            
            smd.put('file', filename)
            smd.put('lemur', lemurname)
            smd.put('date', date)
            smd.put('location', location)
            smd.put('eggs', num_eggs)
            
            self.metadata[filename] = smd
            
        f.close()

