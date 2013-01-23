import sys
import os
from metagenomics.datatypes import Sequence, SampleMetadata, IUPAC

class DataFileError(Exception):
    pass

class DataFile(object) :
    def __init__(self, fname, extension) :
        self.pathname = os.path.dirname(fname)
        self.filename = os.path.basename(fname)
        self.extension = extension

        if not self.__exists() :
            raise DataFileError("'%s' does not exist" % self.get_filename())

    @property
    def name(self) :
        return self.get_filename()

    def __exists(self) :
        return os.path.isfile(self.get_filename())

    def get_filename(self) :
        if self.pathname != "" :
            return self.pathname + os.sep + self.filename
        return self.filename

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

    def __len__(self) :
        return self.__num_states

class ParseError(Exception) :
    pass

class FastqFile(DataFile) :
    SEQID = 0
    SEQ = 1
    QUALID = 2
    QUAL = 3

    def __init__(self, fname) :
        super(FastqFile, self).__init__(fname, ".fastq")
        self._filehandle = None #open(self.get_filename())
        self._state = FastqFile.SEQID
        self._linenum = 0

        self._validators = {
                FastqFile.SEQID  : self.__validate_seqid,
                FastqFile.SEQ    : self.__validate_sequence,
                FastqFile.QUALID : self.__validate_qualid,
                FastqFile.QUAL   : self.__validate_qualities
            }

        self._current = {  
                FastqFile.SEQID  : None,
                FastqFile.SEQ    : None,
                FastqFile.QUALID : None,
                FastqFile.QUAL   : None 
            }

    def __validate_seqid(self, s) :
        if not (s.startswith('@') or s.startswith('>')) :
            raise ParseError("%s : expected line %d to start with a @ or > (started with %s)" % \
                    (self.get_filename(), self._linenum, s[0]))

    def __validate_sequence(self, s) :
        uniq = set(s)

        for i in uniq :
            if i not in IUPAC.codes :
                raise ParseError("%s : line %d contained an invalid UIPAC code (%s)" % \
                        (self.get_filename(), self._linenum, i))

    def __validate_qualid(self, s) :
        if not s.startswith('+') :
            raise ParseError("%s : expected line %d to start with a +" % \
                    (self.get_filename(), self._linenum))

    def __validate_qualities(self, s) :
        uniq = set(s)

        for i in uniq :
            try :
                Sequence.quality_to_int(i)
        
            except ValueError, ve :
                raise ParseError("%s : line %d contained an invalid quality value (%s)" % \
                        (self.get_filename(), self._linenum, i))

    def __validate(self, s, state) :
        self._validators[state](s)

    def __iter__(self) :
        return self

    def next(self) :
        return self.read()

    def open(self) :
        if self._filehandle :
            self._filehandle.close()

        self._filehandle = open(self.get_filename())

    def close(self) :
        if self._filehandle :
            self._filehandle.close()

    def seq(self) :
        duplicates = 1
        if "NumDuplicates" in self._current[FastqFile.SEQID] :
            data = self._current[FastqFile.SEQID].split()
            duplicates = int(data[1].split('=')[1])

        tmp = Sequence(self._current[FastqFile.SEQ], 
                None if self._current[FastqFile.QUAL] == "" else self._current[FastqFile.QUAL])
        
        # hack, maybe make more documented
        tmp.id = self._current[FastqFile.SEQID]
        tmp.duplicates = duplicates

        return tmp

    def read(self) :
        for line in self._filehandle :
            line = line.strip()

            self._linenum += 1

            if line == "" :
                continue

            # both '@' and '>' are legitimate quality scores
            # but '+' is a genuine delimiter
            if line.startswith('+') :
                self._current[FastqFile.QUALID] = line
                self._state = FastqFile.QUAL
                continue

            if self._state == FastqFile.SEQID :
                self.__validate(line, self._state)

                self._current[FastqFile.SEQID] = line
                self._current[FastqFile.SEQ] = ""
                self._current[FastqFile.QUALID] = ""
                self._current[FastqFile.QUAL] = ""

                self._state = FastqFile.SEQ

            elif self._state == FastqFile.SEQ :
                # if we are reading a fasta file
                if line.startswith('>') :
                    tmp = self.seq()
                    self._current[FastqFile.SEQID] = line
                    self._current[FastqFile.SEQ] = ""
                    return tmp

                self._current[FastqFile.SEQ] += line

            elif self._state == FastqFile.QUAL :
                self.__validate(line, self._state)

                self._current[FastqFile.QUAL] += line

                if len(self._current[FastqFile.SEQ]) == len(self._current[FastqFile.QUAL]) :
                    self._state = FastqFile.SEQID
                    return self.seq()

        if self._current[FastqFile.SEQID].startswith('>') :
            tmp = self.seq()
            self._current[FastqFile.SEQID] = ""
            return tmp

        raise StopIteration

class MetadataReader(object) :
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
            
            smd = SampleMetadata()
            
            smd.put('file', filename)
            smd.put('lemur', lemurname)
            smd.put('date', date)
            smd.put('location', location)
            smd.put('eggs', num_eggs)
            
            self.metadata[filename] = smd
            
        f.close()

