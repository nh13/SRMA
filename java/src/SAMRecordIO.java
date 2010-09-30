package srma;

import java.util.*;
import java.io.*;
import srma.AlignRecordComparator;
import net.sf.samtools.*;
import net.sf.samtools.util.*;
import net.sf.picard.sam.*;
import net.sf.picard.io.IoUtil;

public class SAMRecordIO 
{
    // Public
    public SAMFileHeader mergedHeader;
    public SAMProgramRecord programRecord;

    // Private
    private List<SAMFileReader> readers;
    private List<SAMFileHeader> readersHeaders;
    private List<SAMFileWriter> writers;
    private List<CloseableIterator<SAMRecord>> recordsIters = null;
    private List<AlignRecord> buffer = null; // should be one per input file

    public SAMRecordIO(List<File> inputs, List<File> outputs, String programVersion, boolean useRanges, SAMSequenceDictionary referenceDictionary)
        throws Exception
    {
        ListIterator<File> inputsIter = null;
        ListIterator<File> outputsIter = null;
        ListIterator<SAMFileReader> readersIter = null;

        this.readers = new ArrayList<SAMFileReader>();
        this.readersHeaders = new ArrayList<SAMFileHeader>();
        this.writers = new ArrayList<SAMFileWriter>();

        programVersion = new String("srma-" + programVersion); // append "srma-" so we know it was srma

        inputsIter = inputs.listIterator();
        if(1 < outputs.size()) { // to multiple files 
            outputsIter = outputs.listIterator();
            if(outputs.size() != inputs.size()) {
                throw new Exception("There must be the same # of inputs as outputs");
            }
        }
        while(inputsIter.hasNext()) {
            File file = inputsIter.next();
            SAMFileReader fileReader = null;
            SAMFileHeader fileHeader = null;

            IoUtil.assertFileIsReadable(file);

            fileReader = new SAMFileReader(file, true);
            if(useRanges && !fileReader.hasIndex()) {
                throw new Exception("BAM files and BAM indexes when using the RANGE or RANGES option"); 
            }
            fileHeader = fileReader.getFileHeader();
            if(!checkHeaderAgainstReferenceDictionary(fileHeader, referenceDictionary)) {
                throw new Exception("FASTA sequence dictionary and SAM/BAM file dictionary are in different orders");
            }

            this.programRecord = fileHeader.getProgramRecord("srma");

            if(null == programRecord) { // create a new one
                this.programRecord = new SAMProgramRecord("srma");
                this.programRecord.setProgramVersion(programVersion);
                fileHeader.addProgramRecord(this.programRecord);
            }
            else if(0 != programVersion.compareTo(programRecord.getProgramVersion())) { // new version, but srma exists
                this.programRecord = fileHeader.createProgramRecord();
                this.programRecord.setProgramVersion(programVersion);
            }

            // Always set to coordinate sorted
            fileHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);

            this.readers.add(fileReader);
            this.readersHeaders.add(fileHeader);
            if(1 < outputs.size()) { // to multiple files 
                this.writers.add(new SAMFileWriterFactory().makeSAMOrBAMWriter(fileHeader, true, outputsIter.next())); 
            }
        }

        // Merge headers
        final SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(SAMFileHeader.SortOrder.coordinate, this.readersHeaders, true);
        this.mergedHeader = headerMerger.getMergedHeader();
        // Always set to coordinate sorted
        this.mergedHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        if(0 == outputs.size()) { // to STDOUT
            this.writers.add(new SAMFileWriterFactory().makeSAMWriter(this.mergedHeader, true, System.out));
        }
        else if(1 == outputs.size()) { // one output file
            this.writers.add(new SAMFileWriterFactory().makeSAMOrBAMWriter(this.mergedHeader, true, outputs.get(0)));
        }

        // Default iterators
        this.recordsIters = new ArrayList<CloseableIterator<SAMRecord>>();
        readersIter = this.readers.listIterator();
        while(readersIter.hasNext()) {
            this.recordsIters.add(readersIter.next().iterator());
        }

    }

    private boolean checkHeaderAgainstReferenceDictionary(SAMFileHeader header,  SAMSequenceDictionary referenceDictionary)
    {
        int i;
        SAMSequenceDictionary headerDict;
        if(null == header || null == referenceDictionary) {
            return true;
        }

        headerDict = header.getSequenceDictionary();
        if(headerDict.size() != referenceDictionary.size()) {
            return false;
        }
        for(i=0;i<headerDict.size();i++) {
            // check name and length
            SAMSequenceRecord headerRec = headerDict.getSequence(i);
            SAMSequenceRecord referenceRec = referenceDictionary.getSequence(i);

            if(!headerRec.getSequenceName().equals(referenceRec.getSequenceName()) 
                    || headerRec.getSequenceLength() != referenceRec.getSequenceLength()) 
            {
                return false;
            }
        }

        return true;
    }

    private void initBuffer()
    {
        ListIterator<CloseableIterator<SAMRecord>> iter = null;
        int fileIndex = 0;

        this.buffer = new LinkedList<AlignRecord>();

        iter = this.recordsIters.listIterator();
        while(iter.hasNext()) {
            CloseableIterator<SAMRecord> recordsIter = iter.next();
            if(recordsIter.hasNext()) {
                SAMRecord rec = recordsIter.next();
                // search for correct position - implemented for simplicity
                this.addToBuffer(rec, fileIndex);
            }
            fileIndex++;
        }
    }

    private void addToBuffer(SAMRecord rec, int fileIndex) 
    {
        // search for correct position - implemented for simplicity
        // TODO: implement binary search for longer lists
        int i;
        for(i=0;i<this.buffer.size();i++) {
            AlignRecord bufferRec = this.buffer.get(i);
            if(rec.getReferenceIndex() < bufferRec.record.getReferenceIndex() ||
                    (rec.getReferenceIndex() == bufferRec.record.getReferenceIndex() &&
                     rec.getAlignmentStart() <= bufferRec.record.getAlignmentStart())) 
            {
                break;
            }
        }
        this.buffer.add(i, new AlignRecord(rec, null, fileIndex));
    }

    public boolean hasNextAlignRecord()
    {
        if(0 == this.buffer.size()) {
            return false;
        }
        else {
            return true;
        }
    }

    public AlignRecord getNextAlignRecord()
    {
        AlignRecord ar = null;

        if(this.hasNextAlignRecord()) {
            ar = this.buffer.remove(0);
            CloseableIterator<SAMRecord> recordsIter = recordsIters.get(ar.fileIndex);
            if(recordsIter.hasNext()) {
                this.addToBuffer(recordsIter.next(), ar.fileIndex);
            }
        }
        return ar;
    }

    public void query(String sequenceName, int startPosition, int endPosition)
    {
        ListIterator<SAMFileReader> readersIter = null;
        ListIterator<CloseableIterator<SAMRecord>> recordsItersIter = null;

        readersIter = this.readers.listIterator();
        recordsItersIter = this.recordsIters.listIterator();

        // Close all
        while(readersIter.hasNext()) {
            SAMFileReader reader = readersIter.next();
            CloseableIterator<SAMRecord> recordIter = recordsItersIter.next();
            if(reader.hasIndex()) {
                recordIter.close();
                recordsItersIter.set(reader.query(sequenceName, startPosition, endPosition, false));
            }
        }

        this.initBuffer();
    }

    public void output(AlignRecord rec)
    {
        if(1 == this.writers.size()) {
            this.writers.get(0).addAlignment(rec.record);
        }
        else {
            this.writers.get(rec.fileIndex).addAlignment(rec.record);
        }
    }

    public void closeAll()
    {
        int i;

        for(i=0;i<readers.size();i++) {
            readers.get(i).close();
        }
        for(i=0;i<writers.size();i++) {
            writers.get(i).close();
        }
        for(i=0;i<recordsIters.size();i++) {
            recordsIters.get(i).close();
        }

        buffer = null;
    }
}
