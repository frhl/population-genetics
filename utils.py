
## read in hapmap data
def read_snp(file):
    '''Read a snp file into a pandas dataframe'''
    return(pd.read_table(
        file,
        sep='\s+', # columns are separated by whitespace
        # names of the columns
        names=[None, 'chromosome', 'morgans', 'position', 'ref', 'alt'],
        index_col=0
    ))

SNPs = read_snp(path + 'HapMap3.snp') 

def get_chr_range(chromosome):
    '''Returns the range of positions where SNPs for a chromosome are kept'''
    filt = SNPs.query('chromosome=={}'.format(chromosome))
    start = SNPs.index.get_loc(filt.iloc[0].name)
    stop  = SNPs.index.get_loc(filt.iloc[-1].name) + 1
    return(start, stop)

def read_geno(file):
    '''Reads a geno file into a masked numpy matrix'''
    return(np.genfromtxt(
        file,               # the file
        dtype='uint8',      # read the data in as 1-byte integers
        delimiter=1,        # 1-byte width data
        missing_values=9,   # 9 indicates missing data
        usemask=True        # return a masked array
    ))

def read_geno_pop_chr(pop, chromosome):
    '''Reads a slice of a geno file into a masked numpy matrix'''
    f = open(path + pop + '.geno')      # open the file
    (start, stop) = get_chr_range(chromosome)
    s = it.islice(f, start, stop) # slice the file only keeping SNPs of chr
    return read_geno(s)


## functions to calculate stats

def calc_af(geno):
    ''' calculated the allele frequence of a geno file'''
    
    return(geno.mean(axis=1).filled(-1)/2)

def calc_maf(geno):
    ''' calculated mino allele frequence of a geno file'''
    af = calc_af(geno)
    return(1 - af[af > 0.5])
   

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
    
    '''Use hapmap to calculate genetic divergence between population on specifc chromosomes.'''
    
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

def calc_att(M, t):
    
    ''' 
    calculate Cochran–Armitage test for trend within contingency table (M)
    
    M: Contingency table
    t: weights, e.g. [0,1,2]
    
    returns standardized Z-scores
    
    '''
    
    # Get row/col sums
    R = np.sum(M, axis = 1)
    C = np.sum(M, axis = 0)
    N = np.sum(M)
    k = len(t)
    
    # calculate test-statistic
    T = sum([t[i]*(M[0,i]*R[1] - M[1,i]*R[0]) for i in range(0,k)])

    # calc variance by decomposition
    i = 1
    v1 = sum([(t[i]**2)*C[i]*(N-C[i]) for i in range(0,k)])
    v2 = sum([t[i]*t[j]*C[i]*C[j] for i in range(0,k-1) for j in range(i+1,k)])
    q1 = (R[0]*R[1])/N
    q2 = (v1 - 2*v2)
    
    return(T/math.sqrt(q1*q2))




