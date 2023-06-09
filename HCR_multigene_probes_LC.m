% Generates probes for multiplex experiment from multiple genes
% Requires genelist.csv as input with names of genes as they appear in folder
% Needs RNA fasta files to be downloaded from NCBI as: genename_datasets/ncbi_dataset/data/rna.fna
% Change species on line 26 as required - this is for Blast function only 

clear all;

probelength = 52;
spacerlength = 50;

GC_cutoff = [40 65];

genelist = readtable("/Users/cubittl/Documents/Gene_fastas/genelist.csv");  %needs to have Header in cell A1 and then list below 

str = string(genelist{:,1});

for a = 1:length(str)
    fasta_name=strcat("/Users/cubittl/Documents/Gene_fastas/",str{a},"_datasets/ncbi_dataset/data/rna.fna") %this should be the gene fasta file, str{a} is gene name
   
    genename = str{a}
    
    intron="false"
    outfasta=strcat("/Users/cubittl/Documents/Gene_fastas/Output/output_",datestr(now,1), ".fna") %filepath for output file containing probes
    
    finalfasta3=strcat("/Users/cubittl/Documents/Gene_fastas/Final/3probelist_",datestr(now,1),".fna") %filepath for final 3' probe list
    finalfasta5=strcat("/Users/cubittl/Documents/Gene_fastas/Final/5probelist_",datestr(now,1),".fna") %filepath for final 5' probe list
    
    Species='Homo sapiens';

    [Header, Sequence] = fastaread(fasta_name);
    if isa(Sequence, 'char')
        Sequence = {Sequence}
    end

    Sequence = char(Sequence);

    if isa(Sequence, 'char')
        Sequence = string(Sequence); 
    end

    List = {};
    ListTS2 = {};
    probenum = 1;

    sequence_features = 0;
    blast_failure = 0;
    correct = 0;
    total=0;

    for y=1:length(Sequence)
        y
        S = char(Sequence(y));
        S

        tail = probelength + spacerlength + probelength +2; %index to where the tail of the probe design is

        while tail +probelength < length(S)                

            TempSeq = S(tail - probelength - spacerlength: tail-spacerlength);  %define temp sequence to try as probe
            TempSeq2 = S(tail - probelength - spacerlength - probelength - 1: tail - spacerlength + probelength + 1); %define temp seq plus 52 bases either side as flanking regions to avoid in future probes
            total = total + 1;

    %         notSbefore = char(extractBetween(Sequence,1,tail - probelength - spacerlength))     %defines sequence from beginning to start of TempSeq
    %         notSafter = char(extractBetween(Sequence,tail-spacerlength,strlength(Sequence)))    %defines sequence from end of TempSeq to end of fasta

            hits = false;
            blasthits=false;
            
            if any(TempSeq==' ')
                hits = true;
            end
            
            %check for repeats of same base
            if seqwordcount(TempSeq,'CCCCC') > 0
                hits = true;
            end
            if seqwordcount(TempSeq,'GGGGG') > 0
                hits = true;
            end
            if seqwordcount(TempSeq,'AAAAA') > 0
                hits = true;
            end
            if seqwordcount(TempSeq,'TTTTT') > 0
                hits = true;
            end

            %BSAI
            %restriction site exclusion
            if seqwordcount(TempSeq,'GGTCTC') > 0
              hits = true;
            end
            if seqwordcount(TempSeq,'GAGACC') > 0
              hits = true;
            end

            %check GC
            if hits == false
                seqProperties = oligoprop(TempSeq, 'HPBase', 7, 'HPLoop', 6, 'DimerLength', 10);
                if seqProperties.GC < GC_cutoff(1)
                    hits = true;
                end
            end    

            if seqProperties.GC > GC_cutoff(2)
                hits = true;
            end

            %check for hairpins or dimers
            if ~isempty(seqProperties.Hairpins) || ~isempty(seqProperties.Dimers)
                hits = true;
            end

            % check lowercase for repeats! can include 0-5 lowercase bases
            if numel(TempSeq) - nnz(TempSeq == upper(TempSeq)) > 5
                hits = true;
                TempSeq;
            end

            if length(TempSeq) < probelength
                TempSeq
                hits = true;
            end

            %check sequence prior to tempseq for similarity, remove if same sequence present 
    %         if contains(notSbefore, TempSeq)
    %            hits = true
    %         end

            %check sequence after tempseq for similarity, remove if same sequence present
    %         if contains(notSafter, TempSeq)
    %            hits = true
    %         end

            %check whether TempSeq is within flanking regions of other probes
            for j = 1:max(size(ListTS2))
                if contains(ListTS2{j}, TempSeq)
                    hits = true;
                end
            end

            %check whether other probes are within flanking regions of other probes 
            for k = 1:max(size(List))
                if contains(TempSeq2, List{k})
                    hits = true;
                end
            end

% %          %only blast if sequence properties are good
            if hits == false 
%               
                failcounter=0;
                while 1
                    failcounter=failcounter+1;
                    disp('Trying to BLAST sequence.')
                    try
                        [RID1,ROTE] = blastncbi(TempSeq,'blastn','Database','refseq_rna');
                        blastout = getblast(RID1,'WaitTime',ROTE);
                        break
                    catch me
                        if failcounter>25
                            disp('Exceeded maximum tries for BLAST. Sorry!')
                            throw(me)
                        end
                        disp('BLAST failed. Trying again...')
                    end
                end
                disp('BLAST success!')
                for i = 1:length(blastout.Hits)
                for i = 1:numel(blastout.Hits)
                for i = 1:length(blastout.Hits)
                    if ~contains(lower(blastout.Hits(i).Definition), 'predicted') %ignore predicted
                        if contains(lower(blastout.Hits(i).Definition), lower(Species)) %
                            if ~contains(lower(blastout.Hits(i).Definition), lower(genename)) %ignore same gene
                                if isempty(strfind(blastout.Hits(i).Definition, 'RIKEN'))
                                    for j = 1:max(size(blastout.Hits(i).Hsps))
                                        if all(blastout.Hits(i).Hsps(j).Frame==[1,1]) %Look for only Plus/Plus: Frame is [1,1] when plus plus, [1,-1] for plus minus, etc.
                                            if blastout.Hits(i).Hsps(j).Identities >= 25 %25 bases of identity or greater
                                                blastout.Hits(i).Hsps(j).Identities;
                                                if blastout.Hits(i).Hsps(j).QueryIndices(1) < 12 &&  blastout.Hits(i).Hsps(j).QueryIndices(2)>20 %well centered
                                                    alignment = blastout.Hits(i).Hsps(j).Alignment;
                                                    blastout.Hits(i).Definition
                                                    blasthits = true;
                                                    blastout.Hits(i).Hsps(j)
                                                    TempSeq;
                                                end
                                            end
                                        end
                                    end
                                end    
                            end
                        else
                            blastout.Hits(i);
                        end
                    end
                end
                blastres{probenum} = blastout;
                end
                end
            
                
            if blasthits==true       
                tail=tail+8;    
                blast_failure = blast_failure + 1;    
            elseif hits == true
                tail = tail + 5; %shift window by 2 if probes not found
                sequence_features = sequence_features + 1;
            else
                probeindx{probenum} = tail - probelength - spacerlength;
                tail = tail + probelength + spacerlength;
                List{probenum} = TempSeq;
                ListTS2{probenum} = TempSeq2;
    %             ListnotSbefore{probenum} = notSbefore;
    %             ListnotSafter{probenum} = notSafter;
                probenum = probenum + 1;
                TempSeq;
                correct = correct + 1;
            end
            end
        end
      
     
    %reverse comp the TempSeq to form the prove sequence and append to output

    for ii = 1:max(size(List))
        RevList{ii} = seqrcomplement(List{ii});
        fasta.Sequence = RevList{ii};
        fasta.Header = [str{a} '_probe' num2str(ii) ' index ' num2str(probeindx{ii})];
        fastawrite(outfasta, fasta)   
    end
    
   
end    
end

% Next section makes probes into final version with sequences appended
%To change between amplifiers, alter the sequences in lines 235 and 236 as below
    %B1:    probe3 = GAAGAGTCTTCCTTTACG probe5 = gAggAgggCAgCAAACgg
    %B2:    probe3 = ATCATCCAgTAAACCgCC probe5 = CCTCgTAAATCCTCATCA
    %B3:    probe3 = CCACTCAACTTTAACCCg probe5 = gTCCCTgCCTCTATATCT
    %B4:    probe3 = TCTCACCATATTCgCTTC probe5 = CCTCAACCTACCTCCAAC
    %B5:    probe3 = CTACCCTACAAATCCAAT probe5 = CTCACTCCCAATCTCTAT
    
[Header2, Sequence2] = fastaread(outfasta);
for l = 1:length(Sequence2)
               
    n = char(Sequence2(l));
    probe3 = strcat(string(n(1:25)),'TA','TCTCACCATATTCgCTTC'); % 3 prime probe structure Note:this is B1 amplifier sequence
    probe5 = strcat('CCTCAACCTACCTCCAAC','AA',string(n(28:end))); % 5 prime probe structure Note:this is B1 amplifier sequence

    m = char(Header2(l));                               %making header for fasta from initial probe list for final probes
    o = strcat(m, '_3');                                %3'probe header
    p = strcat(m, '_5');                                %5'probe header
    
    fasta.Sequence = probe3;                            
    fasta.Header = [o]
    fastawrite(finalfasta3, fasta);                     %append 3' probe to fasta file

    fasta.Sequence = probe5;                            
    fasta.Header = [p]
    fastawrite(finalfasta5, fasta);                     %append 5' probe to fasta file
        
end 
 
    clear probeindx;
    clear List;

