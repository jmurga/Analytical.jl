# Processing raw data
## muti-FASTA files
We included some tools to process multi-FASTA files into unfolded SFS and divergence. The function [`Analytical.uSfsFromFasta`](@ref) need three files: a reference file to degenerate the sequence, a multi-FASTA file to process the polymorphism and a outgroup sequence to process the divergence.

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

## VCF files
To process VCF files you must include the degeneracy information at the INFO field. It could be easily done through software like SnpEff. To process the divergence we require the reference and the outgroup fasta file as well.
