#./srma -i DH10B.cs.bam -o tmp.bam -r DH10B.fa -R DH10B_WithDup_FinalEdit:1-1000000;
#./srma -i DH10B.cs.bam -o tmp.bam -r DH10B.fa -R DH10B_WithDup_FinalEdit:900000-999994
#./srma -i DH10B.cs.bam -o tmp.bam -r DH10B.fa -R DH10B_WithDup_FinalEdit:1-2000000;
#java -Xmx4G -jar /Users/nilshomer/git/srma/build/jar/srma-0.1.7.jar I=DH10B.cs.bam O=tmp.bam R=DH10B.fa RANGE=DH10B_WithDup_FinalEdit:1-2000000;

#valgrind --leak-check=yes --log-file=log.txt ./srma -i DH10B.cs.bam -o tmp.bam -r DH10B.fa -R DH10B_WithDup_FinalEdit:1-500000;
