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
    private List<SAMFileWriter> writers;
    private List<CloseableIterator<SAMRecord>> recordsIters = null;
    private List<AlignRecord> buffer = null; // should be one per input file
    boolean inputClosed = true;

    public SAMRecordIO(List<File> inputs, List<File> outputs, String programVersion)
        throws Exception
    {
        ListIterator<SAMFileReader> readersIter = null;

        this.init(inputs, outputs, programVersion);

        // Open output iterators ...
        this.recordsIters = new ArrayList<CloseableIterator<SAMRecord>>();
        readersIter = this.readers.listIterator();
        while(readersIter.hasNext()) {
            this.recordsIters.add(readersIter.next().iterator());
        }
        this.initBuffer();
    }

    public SAMRecordIO(List<File> inputs, List<File> outputs, String programVersion, String sequenceName, int startPosition, int endPosition)
        throws Exception
    {
        ListIterator<SAMFileReader> readersIter = null;

        this.init(inputs, outputs, programVersion);

        // Open output iterators ...
        this.recordsIters = new ArrayList<CloseableIterator<SAMRecord>>();
        readersIter = this.readers.listIterator();
        while(readersIter.hasNext()) {
            this.recordsIters.add(readersIter.next().query(sequenceName, startPosition, endPosition, false));
        }
        this.initBuffer();
    }

    private void init(List<File> inputs, List<File> outputs, String programVersion)
        throws Exception
    {
        ListIterator<File> inputsIter = null;
        ListIterator<File> outputsIter = null;

        this.readers = new ArrayList<SAMFileReader>();
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
            fileHeader = fileReader.getFileHeader();
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
            if(1 < outputs.size()) { // to multiple files 
                this.writers.add(new SAMFileWriterFactory().makeSAMOrBAMWriter(fileHeader, true, outputsIter.next())); 
            }
        }

        // Merge headers
        final SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(this.readers, SAMFileHeader.SortOrder.coordinate, true);
        this.mergedHeader = headerMerger.getMergedHeader();
        // Always set to coordinate sorted
        this.mergedHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        if(0 == outputs.size()) { // to STDOUT
            this.writers.add(new SAMFileWriterFactory().makeSAMWriter(this.mergedHeader, true, System.out));
        }
        else if(1 == outputs.size()) { // one output file
            this.writers.add(new SAMFileWriterFactory().makeSAMOrBAMWriter(this.mergedHeader, true, outputs.get(0)));
        }
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
        inputClosed = false;
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

        // Close all
        this.closeInput();

        // Open new output iterators 
        this.recordsIters = new ArrayList<CloseableIterator<SAMRecord>>();
        readersIter = this.readers.listIterator();
        while(readersIter.hasNext()) {
            this.recordsIters.add(readersIter.next().query(sequenceName, startPosition, endPosition, false));
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

    public void closeInput()
    {
        ListIterator<CloseableIterator<SAMRecord>> iter = null;

        // Close all
        iter = this.recordsIters.listIterator();
        while(iter.hasNext()) {
            iter.next().close();
        }

        buffer = null;
        inputClosed = true;
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
        if(!inputClosed) {
            for(i=0;i<recordsIters.size();i++) {
                recordsIters.get(i).close();
            }
        }

        buffer = null;
    }
}
