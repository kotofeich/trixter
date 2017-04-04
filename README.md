TRIXTER is a python utility for managing breakpoints in assemblies.

### INPUT

Input is required as file with synteny blocks in a simple format:<br />
> Block #{number1}<br />
> Specie_id1.Seq_id1_1  Strand1_1  Start2_1   End1_1<br /> 
> Specie_id2.Seq_id2_1  Strand2_1  Start2_1   End2_1<br />  
> Block #{number2}<br />
> Specie_id1.Seq_id1_1  Strand1_1  Start2_1   End1_1<br /> 
> Specie_id2.Seq_id2_1  Strand2_1  Start2_1   End2_1<br />  
> ...<br />

For calculating rearrangements two species must be specified that correspond to Specie_id names.

### TYPES OF ANALYSIS

There are three ways to use TRIXTER

#### Scan the blocks for rearrangements
Here we use the following terminology for rearrangements:<br /> 
1. types of relocation - change of genomic position <br /> 
- translocation of a genomic segment - r. to another chromosome <br /> 
- transposition of a genomic segment - r. to the same chromosome <br /> 
2. inversion - change of strand in a genomic segment <br /> 
3. duplications - copying of a genomic segment <br /> 
4. deletions - removal of a genomic segment <br /> 
We will make the following assumptions: insertion of a new sequence is a pretty seldom event and can be neglected, also 
SNV, CNV, and small indels can be neglected, as we keep a large scale point of view. 
In such a system of operations over the ancestral genome, all the genomes of desendants represented by syntenic blocks 
can be expressed over for this operations.
Here we focus on first two (three) types of rearrangemetns: relocations and inversions. 
All the duplications and insertions are ignored during the analysis.

Disregard of the number of species present in the set of synteny blocks rearrangments can be calculated between a pair 
of species: from one to the other. Thus talking about relocations we report them in the second specie given in the arguments 
with regard to the order of blocks in first specie.

Example:

    breakpoints_analyzer.py --report_translocations —-species Cat Dog --file blocks_coords.txt
   
#### Identify breakpoints among the set of genomes 

All the genomes present in the blocks file are considered.

Example: 

    breakpoints_analyzer.py --classify_breakpoints --file blocks_coords.txt
    
    breakpoints_analyzer.py --classify_breakpoints --print_table --file blocks_coords.txt
    
#### Print out the genome blocks

The blocks file is parsed and the blocks comprising the specified specie chromosomes are printed out.

Example:

    breakpoints_analyzer.py --print_out_genomes Cat --file blocks_coords.txt

### Output




