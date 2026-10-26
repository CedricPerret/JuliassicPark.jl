#***********************************************
#*** Genotype to phenotype mapping
#***********************************************

#--- Trait is average of genotype
"""
    average_mapping(population::AbstractMatrix{<:Real})

Maps a genotype (2D matrix of alleles) to a phenotype by taking the average across each matrix.
"""
average_mapping(genotype::AbstractMatrix{<:Real}) = mean(genotype)
average_mapping(genotypes::Tuple) = map(average_mapping, genotypes)

#--- Trait is sum of genotype
"""
    additive_mapping(population::AbstractMatrix, allele_eff)

Maps a boolean genotype matrix to a phenotype by computing the additive effect `allele_eff × (2A - 2n)`,
where `A` is the number of true alleles and `n` is the number of loci.
"""
additive_mapping(genotype::AbstractMatrix{Bool}, allele_eff) = allele_eff * (2 * sum(genotype) - 2 * size(genotype, 1))
additive_mapping(genotypes::Tuple, allele_eff) = map(g -> additive_mapping(g, allele_eff), genotypes) 
