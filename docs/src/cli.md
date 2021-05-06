# Command line interface

We develop a Command Line Interfence in case you want to avoid Julia interpreter. You can easly download [abcmk_cli.jl](https://raw.githubusercontent.com/jmurga/Analytical.jl/master/scripts/abcmk_cli.jl). The CLI have different functions to perform the whole pipeline as explained in [Infering the rate and strength of adaptation](@ref) section.

```bash
julia abcmk_cli.jl  
```

```bash
See --help of each command for usages  
  rates  
  parseTgpData  
  parseDgnData  
  summaries  
  abcInference  
  plotMap  
```

To reproduce the examples you can easly:

1. Estimate rates
```bash
time julia abcmk_cli.jl rates --samples 661 --gamNeg "-2000 -200" --gL "1 10" --gH "200 2000" --rho 0.001 --theta 0.001 --solutions 100000 --output /home/jmurga/rates.jld2 --dac 1,2,4,5,10,20,50,100,200,400,500,661,925,1000 --nthreads 22
```

2. Parse data into new folder
```bash
julia abcmk_cli.jl parseData --analysisFolder /home/jmurga/tgpData
```

3. Estimate summary statistics
```bash
julia abcmk_cli.jl summaries --analysisFolder /home/jmurga/tgpData
```

4. Estimate summary statistics
```bash
julia abcmk_cli.jl summaries --analysisFolder /home/jmurga/tgpData/ --rates /home/jmurga/rates.jld2 --samples 661 --replicas 100 --summstatSize 100000 --dac 2,4,5,10,20,200,661,925 --nthreads 22
```

5. ABC inference
```bash
julia abcmk_cli.jl abcInference --analysisFolder /home/jmurga/tgpData/ --replicas 100 --P 5 --S 9 --tol 0.001 --ABCreg /home/jmurga/ABCreg/src/reg --parallel true --nthreads 22
```

6. ABC inference
```bash
julia abcmk_cli.jl abcInference --analysisFolder /home/jmurga/tgpData/ --replicas 100 --P 5 --S 9 --tol 0.001 --ABCreg /home/jmurga/ABCreg/src/reg --parallel true --nthreads 22
```