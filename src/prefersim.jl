
const T = gsl_rng_default;
const r = rng_alloc(T);

@with_kw struct recipe
	epochs::Array{Int64,1} = [1000]
	N::Array{Int64,1} = [8000]

	θ::Float64=2000
	h::Float64=0.5
	
	s₋::Float64=-457.0
	s₊::Float64=500.0

	dfe::String="point"
	param_one::Float64=1.0
	param_two::Float64=1.0
	s::Array{Float64,1}=-[0.0]
	s_mult::Array{Float64,1}=[1.0]
	prob::Array{Float64,1}=[0.0]

	n_anc::Int64=1000
	burnin_period::Bool=false
	
	relax::Bool=false
	epoch_relaxation::Array{Bool}=fill(false,length(epochs))
	s_relaxation::Float64=0.0
	s_relaxation_threshold::Float64=0.0

	F::Array{Float64,1}=zeros(length(N))
	seed::Int64=rand(1:10^8)

	trajectories::Array{Int64,1} = Int64[]
	trajectories_output::String = ""

	@assert length(N)==length(epochs)  "N and epochs must be equal in length";
	@assert length(s)==length(prob)    "s and probs must be equal in length";
end

@with_kw mutable struct prf_output
	fixations::Int64 = 0;
	loss::Int64=0;
	total_mut::Int64 = 0;
	count_mut::Int64=0;
end

@with_kw mutable struct mutation
	# The frequency of the mutation
	frequency::Float64=0;
	# Count of the mutation
	count_samp::Int64=1; 
	# Selection coefficient
	s::Float64=0.0; 
	# Dominance factor
	h::Float64=0.5; 
	# The generation that mutation has arise
	age::Int64=1;
	# Number of mutation. Tag to avoid similar nodes
	num::Int64=1;
	# Mutation type. 1=neutral; 2=deleterious; 3=strong advantegeous; 4=weakly advantegeous
	type::Int8=0
end

function relax_selection!(mutation_list::LinkedList{mutation}, new_s, sel_threshold::Float64, relax_type::Int64)

	if(!isempty(mutation_list))

		# Index to delete node
		node=mutation_list.node.next;
		while (node.next != mutation_list.node.next)
			item        = node.data;

			freq            = item.frequency;
			s               = item.s;
			h               = item.h;
			freq            = freq_inbreed(s,freq,h);
			num       = item.num;

			if ((s/2) <= sel_threshold && relax_tye == 0)
				item.s = new_s * 2;
			elseif ((s/2) <= sel_threshold && relax_tye == 1)
				item.s = s * new_s * 2;
			end
			node  = node.next
		end
	end
end

function freq_inbreed(sel::Float64,freq::Float64,h::Float64,F::Float64) 
	#The fitnesses of the genotypes are A1A1 = 1; A1A2 = 1-h*(2*s); A2A2 = 1-(2*s), where A2 is the derived allele. Under any conditions that the program is run, the value of s for any segregating site cannot exceed 0.5.
	# parameters struct assume positive values of s to concieve positive selection. Changing to negative value
	s::Float64=-sel
	return (((1.0 - s) * (freq * freq + F * freq * (1.0 - freq))) + ((1.0 - h * s) * freq * (1.0 - freq) * (1.0 - F))) / (((1.0 - freq) * (1.0 - freq) + (1.0 - freq) * F * freq) + ((1.0 - h * s) * 2.0 * freq * (1.0 - freq) * (1.0 - F)) + ((1.0 - s) * (freq * freq + F * freq * (1.0 - freq))))
end

function add_mutation!(mutation_list::LinkedList{mutation},N::Float64,h::Float64,s::Float64,θ::Float64,freq::Float64,dfe::String,param_one::Float64,param_two::Float64,s_mult::Float64,n_anc::Int64,age::Int64,trajectories::Array{Int64,1},relax::Bool=false) 

	num_mut::Float64 = ran_poisson(r,θ/2.0)
	count_mut::Int64 = length(mutation_list)

	for x::Int64=1:num_mut

		if dfe == "point"
			s = s;
		elseif dfe == "gamma"
			##/Use this for gamma distribution in Boyko 2008:
			gamma = ran_gamma(r, param_one, param_two * s_mult);
			# note, this is scale for a Nanc=1000, using the boyko params
			s = - gamma / (n_anc * 2);
		end

        count_mut += 1;
        mut        = mutation(frequency = freq, h = h, s = s*2, count_samp = 0.0, age = age, num = count_mut, type=1)
		push!(mutation_list,mut);
		
		if(!isempty(trajectories) && count_mut in trajectories)
			i = findfirst(isequal(count_mut),trajectories);
			@printf(io, "%d\t%lf\n", trajectories[i], (1.0/N));
		end
	end
end

function add_mutation!(mutation_list::LinkedList{mutation},N::Float64,h::Float64,s::Float64,θ::Float64,freq::Float64,age::Float64) 

	num_mut::Float64 = ran_poisson(r,θ/2.0)
	count_mut::Int64 = length(mutation_list)

	for x::Int64=1:num_mut

		s = s;

        count_mut += 1;
        mut        = mutation(frequency = freq, h = h, s = s*2, count_samp = 0.0, age = age, num = count_mut, type=1)
		push!(mutation_list,mut);

	end
end

function drift_sel!(mutation_list::LinkedList{mutation},N::Float64,F::Float64,trajectories::Array{Int64,1})

	l = 0;
	f = 0;

	if(!isempty(mutation_list))

		# Index to delete node
		node=mutation_list.node.next;
		while (node.next != mutation_list.node.next)
			item        = node.data;

            freq            = item.frequency;
            s::Float64      = item.s;
            h::Float64      = item.h;
            freq::Float64   = freq_inbreed(s,freq,h,F);
            
            num       = item.num;

            #count::Int64    = rand(Binomial(N,freq));
            count::Int64    = ran_binomial(r,freq,N)
            item.count_samp = count;
            freq            = count / N;
            item.frequency  = freq;

			if (freq > 0.0 && freq < 1.0)	
				if(!isempty(trajectories) && num in trajectories)
					i = findfirst(isequal(num),trajectories);
					@printf(io, "%d\t%lf\n", trajectories[i], freq);
				end
				# Move node
				node = node.next
			else
				if(freq == 0)
					l +=1;
				elseif(freq == 1)
					f +=1;
				end
				# Delete node
				deleteat!(mutation_list,node);
				# Move node
				node = node.next
			end
		end
		return(l,f)
	else
		return(l,f)
	end
end

function burnin(param::recipe)

	@unpack N,θ,s,h,dfe,param_one,param_two,s_mult,n_anc,trajectories,relax = param;

	n_size = N[1]-1;
	s_size = N[1]-1;
	sel    = s[1];

	# theoretical_sfs     = zeros(n_size,2);
	number_of_mutations = 0;
	
	mutation_list = LinkedList{mutation}()

	f(x,j) = f_of_q_lambda(x,j,n_size,s_size,sel,θ*2)
	@inbounds for j::Int64=2:n_size
		m = LinkedList{mutation}()
		
		res,err = quadgk(x->f(x,j),0,1)
		number_of_mutations += res;
		age = n_size * 10;
		
		add_mutation!(m,Float64(n_size),h,sel,res,j/n_size,dfe,param_one,param_two,s_mult[1],n_anc,age,trajectories,relax);
		append!(mutation_list,m);
	end
	return(mutation_list);
end

function f_of_q_lambda(x::Float64,j::Int64,N_burnin::Int64,sample_size::Int64,point_sel::Float64,theta::Float64)

	gamma = N_burnin * point_sel;

	if (abs(gamma) > 1.0e-7)
		
		sfs_function_term_one = (1-exp(-2*N_burnin*(-point_sel)*(1-x)))/(1-exp(-2*N_burnin*(-point_sel)));
		sfs_function_term_two::Float64 = (2/(x*(1-x)));

		binomial_value::Float64 = ran_binomial_pdf(j,x,sample_size);

		f = theta/2 * sfs_function_term_one * sfs_function_term_two * binomial_value;

	else
		f = theta/2 * 2/x * ran_binomial_pdf(j,x,sample_size);
	end
	return  f;
end

function sfs(mutation_list::LinkedList{mutation},sample_size)

	if(!isempty(mutation_list))

		mut  = zeros(length(mutation_list)); 
		#=neut = zeros(length(mutation_list));
		del  = zeros(length(mutation_list));
		adv  = zeros(length(mutation_list));=#
		
		i    = 1;
		# Index to delete node
		node = mutation_list.node.next;

		while (node.next != mutation_list.node.next)
			item        = node.data;

			freq        = item.frequency;
			s           = item.s;

			h           = item.h;
			#samp_count = rand(Binomial(sample_size,freq));
			samp_count = ran_binomial(r,freq,sample_size);
			samp_freq = samp_count/sample_size;
			if (samp_freq > 0 && samp_freq < 1.0)
				mut[i] = samp_freq
				#=if(item.type == 1)
					neut[i] = samp_freq
				elseif(item.type == 2)
					del[i] = samp_freq
				elseif(item.type == 3)
					adv[i] = samp_freq
				end=#
			end
			node = node.next;
			i+=1
		end

		mut     = mut[mut .!= 0]

		b        = collect(1/sample_size:1/sample_size:1)
		v(x,b=b) = searchsortedfirst.(Ref(b), x)
		n    = v(mut)
		
		out      = zeros((sample_size-1,2))
		out[:,1] = collect(1:(sample_size-1))

		@inbounds for i::Int64 in out[:,1]
			out[i,2] = length(n[n .== i])
		end

		out[:,1] = round.(out[:,1]/sample_size,digits=3);

		return out
	else
		out = zeros((sample_size-1,3));
		out[:,1] = round.(out[:,1]/sample_size,digits=3);
		return out
	end
end

function simulate(param::recipe, sample_size::Int64)
	

	@unpack epochs,N,θ,h,s₋,s₊,dfe,param_one,param_two,s,s_mult,prob,n_anc,burnin_period,relax,epoch_relaxation,s_relaxation,s_relaxation_threshold,F,trajectories,trajectories_output,seed = param;

	if !isempty(trajectories)
		@assert length(unique(trajectories .> 0)) == 1  "ID index must be greater than 0";
		io = open(trajectories_output,"w+")
	end
	
	# set seed before start
	rng_set(r, seed)

	epochs = SVector{length(epochs)}(epochs)
	s_mult = SVector{length(s_mult)}(s_mult)
	s = SVector{length(s)}(s)
	F = SVector{length(F)}(F)

	events = length(epochs);
	
	if(length(s)==1)
		s      = s[1];
		s_mult = s_mult[1]
	else
		s = s;
		s_mult = s_mult;
	end
	

	@printf "Demographic History (%i epochs)\n\n" events;
	for e=1:events
		@printf "Ne = %i\tgenerations = %i\tF = %lf\n" N[e] epochs[e] param.F[e];
	end

	mutation_list = LinkedList{mutation}();

	if !burnin_period
		# mutation_list = burnin(mutation_list,param);
		mutation_list = burnin(param);
	else
		mutation_list = LinkedList{mutation}();
	end;
	
	l = 0;
	f = 0;
	N_F_1=0;
	@printf "\n"
	@inbounds for e=1:events
		# Inbreeding Ne

		N_F = N[e] / (1.0 + F[e]); 
		@printf "Currently in epoch = %i ; Mutations before epoch's beginning = %i\n"  e length(mutation_list);

		if epoch_relaxation[e]
			printf("Relaxation in Epoch %i\n", e);
			relax_selection!(mutation_list, s_relaxation, s_relaxation_threshold, relaxation_type);
		end

		@inbounds for g=1:epochs[e]
			# Mutation age
			if e == 1 
				age = g;
			else
				# Add age to previous epochs generations
				age = g + epochs[e];
			end
			x, y    = drift_sel!(mutation_list,N_F,F[e],trajectories);
			add_mutation!(mutation_list,N_F,h,s,θ,1.0/N[e],dfe,param_one,param_two,s_mult,n_anc,age,trajectories,relax);
			l = l + x;
			f = f + y;
		end

		if (e < events) 
			N_F_1 = N[e + 1] / (1.0 + param.F[e + 1]);
			θ = θ * N_F_1 / N_F;
		end
	end

	if @isdefined io
		close(io)
	end

	out = sfs(mutation_list,sample_size)

	return(out,hcat(l,f))
end

function simulate_batch(param::recipe,sample_size::Int64,replicas::Int64)

	mutation_list = LinkedList{mutation}();

	n_params = [deepcopy(param) for i::Int64=1:replicas];
	n_mutations_list = [deepcopy(mutation_list) for i::Int64=1:replicas];
	n_sample_size = [sample_size for i::Int64=1:replicas]


	mut = pmapbatch(simulate,n_params,n_sample_size);
	s = sum(x->x[1][:,2],mut);
	f = sum(x->x[2],mut);

	s = hcat(round.(collect(1:(sample_size-1))/sample_size,digits=3),s);

	return s,f
end

function two_epochs(ne_prior::Union{Int64,Nothing},θ_prior::Union{Float64,Nothing},sample_size::Int64,replicas::Int64;expansion=false,bottleneck=false)

	if(isnothing(θ_prior))
		θ_prior = rand(50:1000,replicas)
	else
		θ_prior = fill(θ_prior,replicas)
	end

	if(isnothing(ne_prior))
		ne_prior = rand(1000:10000,replicas)
	else
		ne_prior = fill(ne_prior,replicas)
	end


	n_sample_size = fill(sample_size,replicas)
	n_params = recipe[];
	
	if expansion
		ne_rand = rand(2:50,replicas)
	elseif bottleneck
		ne_rand = rand(0.01:0.01:1,replicas)
	else
		ne_rand = rand(0.01:0.01:50,replicas);
	end

	ne_two      =  Int64.(trunc.(ne_rand.*ne_prior,digits=0));
	epoch_two   = rand(100:2000,replicas)

	for i::Int64=1:replicas
		push!(n_params,recipe(N=[ne_prior[i],ne_two[i]],epochs=[1,epoch_two[i]],θ=θ_prior[i],burnin_period=false))
	end
	
	mut = pmapbatch(simulate,n_params,n_sample_size);
	
	f(x) = map(y->y[1][:,2],x)
	ps   = hcat(f(mut)...)
	
	return(ps,hcat(ne_prior,ne_two,epoch_two))
end

function three_epochs(ne_prior::Union{Int64,Nothing},θ_prior::Union{Float64,Nothing},sample_size::Int64,replicas::Int64;expansion=false,bottleneck=false)
	
	if(isnothing(θ_prior))
		θ_prior = rand(50:1000,replicas)
	else
		θ_prior = fill(θ_prior,replicas)
	end

	if(isnothing(ne_prior))
		ne_prior = rand(1000:10000,replicas)
	else
		ne_prior = fill(ne_prior,replicas)
	end

	if expansion
		ne_rand       = rand(2:50,replicas)
		ne_rand_three = rand(2:50,replicas)
	elseif bottleneck
		ne_rand       = rand(0.01:0.01:1,replicas)
		ne_rand_three = rand(0.01:0.01:1,replicas)
	else
		ne_rand       = rand(0.01:0.01:50,replicas);
		ne_rand_three = rand(0.01:0.01:50,replicas);
	end

    ne_two        = Int64.(trunc.(ne_rand.*ne_prior,digits=0));
    ne_three      = Int64.(trunc.(ne_rand_three.*ne_two,digits=0));

    epoch_two     = rand(100:1000,replicas)
    epoch_three   = rand(100:1000,replicas)

    n_sample_size = fill(sample_size,replicas)
    n_params      = recipe[];
	for i::Int64=1:replicas
		push!(n_params,recipe(N=[ne_prior[i],ne_two[i],ne_three[i]],epochs=[ne_prior[i]*10,epoch_two[i],epoch_three[i]],θ=θ_prior[i]))
	end

	mut = pmapbatch(simulate,n_params,n_sample_size);

	f(x) = map(y->y[1][:,2],x)
	ps   = hcat(f(mut)...)
	
	return(ps)
end

function exponential(ne_prior::Union{Int64,Nothing},θ_prior::Union{Float64,Nothing},rate::Union{Float64,Nothing},sample_size::Int64,replicas::Int64;expansion=false,bottleneck=false)

	if(isnothing(θ_prior))
		θ_prior = rand(50:1000,replicas)
	else
		θ_prior = fill(θ_prior,replicas)
	end

	if(isnothing(ne_prior))
		ne_prior = rand(1000:10000,replicas)
	else
		ne_prior = fill(ne_prior,replicas)
	end

	if(isnothing(rate))
		rate = rand(0.00005:0.0001:0.02,replicas)
	else
		rate = fill(rate,replicas)
	end

	n_sample_size = fill(sample_size,replicas)
	n_params = recipe[];
	
	if expansion
		ne_rand = rand(2:10,replicas)
	elseif bottleneck
		ne_rand = rand(0.01:0.01:1,replicas)
	else
		ne_rand = rand(0.01:0.01:50,replicas);
	end

	ne_two      =  trunc.(Int64,ne_rand.*ne_prior);
	epoch_two   = rand(100:10000,replicas)
	epoch_exp   = rand(100:1000,replicas)
	
	g(p::Int64,r::Float64,t::Int64) = trunc(Int64,p*(1 + r)^t)
	for i::Int64=1:replicas
		tmp_ne    = g.(ne_two[i],rate[i],1:epoch_exp[i])
		tmp_epoch = fill(1,epoch_exp[i])

		push!(n_params,recipe(N=vcat(ne_prior[i],ne_two[i],tmp_ne),epochs=vcat(1,epoch_two[i],tmp_epoch),θ=θ_prior[i],burnin_period=false))
	end
	
	mut = pmapbatch(simulate,n_params,n_sample_size);
	
	f(x) = map(y->y[1][:,2],x)
	ps   = hcat(f(mut)...)
	
	return(ps,hcat(ne_prior,ne_two,epoch_two,epoch_exp,rate))
end
