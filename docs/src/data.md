# Parsing genomic data
The module includes functions to parse TGP from Uricchio et al. (2019) and DGN from Murga-Moreno et al. (2019). In addition, the module have a function to parse SFS and divergence from multi-FASTA data following Murga-Moreno et al. (2019)

Please to parse raw data into SFS and divergence counts, first download raw files deposited in our repository:  

 - [TGP](https://raw.githubusercontent.com/jmurga/Analytical.jl/master/data/tgp.txt)
 - [DGN Zambia population](https://raw.githubusercontent.com/jmurga/Analytical.jl/master/data/dgnRal.txt)  
 - [DGN Raleigh population](https://raw.githubusercontent.com/jmurga/Analytical.jl/master/data/dgnZi.txt)  

```bash
curl -o analysis/tgp.txt https://raw.githubusercontent.com/jmurga/Analytical.jl/master/data/tgp.txt
```
## Parsing TGP and DGN data manually
Once you have downloaded the files, you can use the function ```Analytical.parseSfs``` to convert the data into SFS and divergence counts. Please check [`Analytical.parseSfs`](@ref) to get more info o execute:

```julia
alpha, sfs, divergence = Analytical.parseSfs(sampleSize = 661, data = "analysis/tgp.txt")
```

To save the data, you can use CSV and DataFrames packages

```julia
using CSV, DataFrames
CSV.write("analysis/tgpSfs.tsv",DataFrame(sfs),delim='\t',header=false))
CSV.write("analysis/tgpDiv.tsv",DataFrame(permutedims(divergence)),delim='\t',header=false))
```

It is possible to directly subset genes IDs using Ensembl or Flybase id. Use a variable of type ```Matrix{String}``` into the argument *geneList*

```julia
ensemblList = CSV.read("analysis/ensemblList.txt",header=false,DataFrame) |> Array

alpha, sfs, divergence = Analytical.parseSfs(sampleSize = 661, data = "analysis/tgp.txt",geneList = ensemblList)
```

If you are going to parse DGN, you need to change the value of the argument *isoline* to *true*. Following the Murga-Moreno et al. (2019) sample size for each population is:

 - Zambia population: 154
 - RAL population: 160

```julia
alpha, sfs, divergence = Analytical.parseSfs(sampleSize = 160, data = "analysis/dgnRal.txt",isolines=true)
```

## Processing muti-FASTA files
We included some tools to process multi-FASTA files into unfolded SFS and divergence. The function [`Analytical.uSfsFromFasta`](@ref) needs three files: a reference file to degenerate the sequence and a multi-FASTA file to process the polymorphism, and an outgroup sequence to process the divergence. Consider downloading DGN data from [John Pool lab](https://www.johnpool.net/) or [PopFly](https://popfly.uab.cat) to obtain these files.

```julia
<!-- sfs, div = Analytical.uSfsFromFasta(
                        file = "/home/jmurga/Downloads/example.fa",
                        reference  = "/home/jmurga/Downloads/ref.fa",
                        outgroup   = "/home/jmurga/Downloads/outgroups.fa",
                        samples    = 160,
                        bins       = 20,
                        codonTable = "standard"
                     ) -->
```

The scripts transforms a multi-FASTA file into a matrix which eliminates the monomorphic sites and process the divergent and polymorphic sites. Only 0-fold and 4-fold sites are analyzed as proxy of neutral and selected alleles. The following matrix correspond to the matrix processed at the above example.

```julia
<!-- 162×153456 Array{Char,2}:
 '0'  '0'  '0'  '0'  '0'  '0'  '4'  …  '0'  '0'  '4'  '0'  '0'  '0'
 'A'  'A'  'G'  'A'  'C'  'C'  'T'     'G'  'G'  'C'  'A'  'G'  'T'
 'A'  'A'  'G'  'A'  'C'  'C'  'T'     'G'  'G'  'C'  'A'  'G'  'T'
 'A'  'A'  'G'  'A'  'C'  'C'  'T'     'G'  'G'  'C'  'A'  'G'  'T'
 'A'  'A'  'G'  'A'  'C'  'C'  'T'     'G'  'G'  'C'  'A'  'G'  'T'
 'A'  'A'  'G'  'A'  'C'  'C'  'T'  …  'G'  'G'  'C'  'A'  'G'  'T'
 'A'  'A'  'G'  'A'  'C'  'C'  'T'     'G'  'G'  'C'  'A'  'G'  'T'
 'A'  'A'  'G'  'A'  'C'  'C'  'T'     'G'  'G'  'C'  'A'  'G'  'T'
 ⋮                        ⋮         ⋱  ⋮                        ⋮
 'A'  'A'  'G'  'A'  'C'  'C'  'T'     'G'  'G'  'C'  'A'  'G'  'T'
 'A'  'A'  'G'  'A'  'C'  'C'  'T'  …  'G'  'G'  'C'  'A'  'G'  'T'
 'A'  'A'  'G'  'A'  'C'  'C'  'T'     'G'  'G'  'C'  'A'  'G'  'T'
 'A'  'A'  'G'  'A'  'C'  'C'  'T'     'G'  'G'  'C'  'A'  'G'  'T'
 'A'  'A'  'G'  'A'  'C'  'C'  'T'     'G'  'G'  'C'  'A'  'G'  'T'
 'A'  'A'  'G'  'A'  'C'  'C'  'T'     'G'  'G'  'C'  'A'  'G'  'T'
 'A'  'A'  'G'  'A'  'C'  'C'  'T'  …  'G'  'G'  'C'  'A'  'G'  'T'
 'A'  'A'  'G'  'C'  'C'  'A'  'T'     'G'  'G'  'C'  'A'  'G'  'T' -->
```

Each column correspond to a sample position. Only fixed divergence and biallelic sites are process. The matrix is iterated through columns to discard the rest of the cases. The matrix is iterated to check the columns independently. In this way we check fixed divergence or biallelic polymorphism. For example at the following array all polymorphic positions (from row 2 to row 161) correspond to 'A' nucleotides where as the outgroup (last row) correspond to 'G'. This columns is addedas to neutral or selected divergence depedending on the first row value:

```julia
<!-- 162×1 Array{Char,1}:
 '0'
 'A'
 'A'
 'A'
 'A'
 'A'
 'A'
 'A'
 ⋮  
 'A'
 'A'
 'A'
 'A'
 'A'
 'A'
 'A'
 'C' -->
```

The next example represent one polymorphic site taking into account the derived allele. At this array, the first 3 polymorphic positions correspond to 'G' where as the rest of polymorphic sites correspond to C. Since 'C' is the ancestral nucleotide (last row) the derived allele frequency will be 3/160:
```julia
<!-- 162×1 Array{Char,1}:
 '4'
 'G'
 'G'
 'G'
 'C'
 'C'
 'C'
 'C'
 ⋮  
 'C'
 'C'
 'C'
 'C'
 'C'
 'C'
 'C'
 'C' -->
```


