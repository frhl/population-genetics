
def calc_fst1(p1, p2, N1, N2):
    
    '''
    Use allele frequencies to calculate genetic diverence
    
            sum( (p1-p2)**2 - ( 1/(2*N1 ) +  1/(2*N2))*p*(1-p))
     fst = -----------------------------------------------------
                          sum(2*p*(1-p))
    
    '''
    p = (p1+p2)/2
    numerator = sum((p1-p2)**2 - (1/(2*N1) + 1/(2*N2))*p*(1-p))
    denominator = sum(2*p*(1-p))
    fst1 = numerator/denominator
    return(fst1)
    


def wrapper_calc_pop_fst1(pop1, pop2, chrom):
    
    '''
    Use hapmap to calculate genetic divergence between population on specifc chromosomes.
    
    '''
    
    # get genotypes
    pop1_geno = read_geno_pop_chr(pop1, chrom)    
    pop2_geno = read_geno_pop_chr(pop2, chrom)
    
    # get sample sizes and allele frequencies
    N1 = pop1_geno.shape[0]
    N2 = pop2_geno.shape[0]
    p1 = calc_af(pop1_geno)
    p2 = calc_af(pop2_geno)
    
    # caluate fst
    fst = calc_fst1(p1, p2, N1, N2)
    return(fst)
    
def is_val(x):
    '''
    checks if is a non-none value
    '''
    return(x != None)

def drop_none(x, y):
    '''
    Drop SNP alleles that are not valid (e.g. None)
    
    '''
    x = np.asarray(x.tolist()) 
    y = np.asarray(y.tolist())
    v = np.vectorize(is_val)(x) & np.vectorize(is_val)(y)
    return(x[v], y[v])

def drop_monomorphic(geno):
    ''' 
    drop monomorphic SNPs (MAF = 0%) from unphased diploid genotype data. 
    
    '''
    af = calc_af(geno)
    return(geno[af != 1])


def r2(a, b):
    '''
    r2 = (E(gA*gA) - E(gA) * E(gB))**2 / var(gA) * var(gB)
    
    calculate LD between a pair of a SNPs.
    
    '''
    a, b = drop_none(a, b)
    num = np.power(np.mean(a * b) - np.mean(a) * np.mean(b), 2)
    denom = np.var(a) * np.var(b)
    return num / denom



