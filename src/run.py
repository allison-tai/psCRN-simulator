import sys, os, configparser, ast, copy
from pqdict import PQDict
import numpy as np
from numpy import array, zeros, log, seterr
from numpy.random import rand, randint
from collections import Counter
from datetime import datetime

# initialize global variables
species_list = []
poly_dict = {}

class Reaction(object):
    def __init__(self, count, r_idx, p_idx, pl_idx, rxn_cat):
        self.count = count
        self.r_idx = r_idx
        self.p_idx = p_idx
        self.pl_idx = pl_idx
        self.rxn_cat = rxn_cat

    def propensity(self, config, volume):
        c = config[self.r_idx]
        if (self.r_idx[0] == self.r_idx[1]): # bi-homogenous
            if (c[0] > 1):
                return c[0]*(c[0]-1)/(2*volume)
            else:
                return c[0]*(0)/(2*volume)
        else: #bi-different
            return (c[0]*c[1])/(volume)
    
    def stoich(self):
        n = len(species_list)
        stoich = zeros(n, dtype=int)
        for sp_idx in self.count:
            stoich[sp_idx] = self.count[sp_idx]
        return stoich

def insert_sp(sp): # returns the index
    try:
        return species_list.index(sp)
    except ValueError:
        species_list.append(sp)
        return len(species_list) - 1

def insert_pl(sp_idx): # returns the index
    pl_idx = poly_dict.get(sp_idx, len(poly_dict))
    if (pl_idx == len(poly_dict)):
        poly_dict[sp_idx] = pl_idx
    return pl_idx

def get_rxncat(reactants, products): #categories: 0 (normal) 1 (push) 2 (pop) 3 (empty) 4 (destroy polymer) 5 (create polymer) 
    rpoly = False
    ppoly = False
    push = False
    pop = False
    for r in reactants:
        if 'bot' in r:
            rpoly = True
        if 'A_' in r:
            pop = True
    for p in products:
        if 'bot' in p:
            ppoly = True
        if 'A_' in p:
            push = True
    if rpoly and ppoly and push:
        return 1 # push
    elif rpoly and ppoly and pop:
        return 2 # pop
    elif rpoly and ppoly:
        return 3 # empty
    elif rpoly:
        return 4 # destroy
    elif ppoly:
        return 5 #create
    else:
        return 0 #normal
    
def new_reaction(reactants_str, products_str): 

    reactants = reactants_str.split(' + ')
    products  = products_str.split(' + ')
    rxn_cat = get_rxncat(reactants, products)

    r_idx = []
    p_idx = []
    pl_idx = -1

    for sp in set(reactants + products):
        sp_idx = insert_sp(sp)
        if ('bot' in sp): #polymer reaction
            pl_idx = insert_pl(sp_idx) # as of now only one polymer per reaction possible
        if sp in reactants:
            r_idx.append(sp_idx)
        if sp in products:
            p_idx.append(sp_idx)
    
    count = Counter(p_idx)
    count.subtract(Counter(r_idx))

    return Reaction(count, r_idx, p_idx, pl_idx, rxn_cat)

def extract_rxns(file):
    protocol_f = open(file, 'r')

    reactions = []
    for line in protocol_f:
        reactions.append(new_reaction(*line.strip().split(' -> ')))

    return reactions

def get_dependencies(reactions):
    dependencies = {}
    for rxn in reactions:
        # dependency list
        dep_list = []
        agent_idx = set(i for i in rxn.count if rxn.count.get(i) != 0)
        # get indices of the reactants/products
        agents = set(species_list[i] for i in agent_idx)
        for rxn2 in reactions:
            # only need indices of reactants for second reaction
            agents2 = set(species_list[i] for i in rxn2.r_idx)
            # small hack: if following our compiler, reactions that only share blanks are never dependent
            if bool(agents & agents2) and agents2 & agents != set('B'):
            # if bool(agents & agents2) # safer version
                dep_list.append(rxn2)
        # add key-value pair (reaction and its dependent reactions) to dictionary
        dependencies[rxn] = dep_list
    return dependencies

def compute_volume(c_0, c_poly):
    pl_total = 0
    for pl in c_poly:
        pl_total = pl_total + sum(pl)
    return sum(c_0) + pl_total

def pl_cases(rxn, react, c_poly): #processes effects of polymer lists of polymer reactions
    if (rxn.rxn_cat <= 4):
    # pick a polymer
        cpl_idx = randint(0, len(c_poly[rxn.pl_idx]))
        # push reaction
        if (rxn.rxn_cat == 1):
            # push
            c_poly[rxn.pl_idx][cpl_idx] += 1
        # pop reaction
        elif (rxn.rxn_cat == 2):
            # if polymer is not empty
            if c_poly[rxn.pl_idx][cpl_idx] > 0:
                # pop
                c_poly[rxn.pl_idx][cpl_idx] -= 1
            # else do nothing
            else:
                react = 0
        # empty reaction
        elif (rxn.rxn_cat == 3):
            if c_poly[rxn.pl_idx][cpl_idx] != 0:
                react = 0
        # destroy reaction
        elif (rxn.rxn_cat == 4):
            if c_poly[rxn.pl_idx][cpl_idx] == 0:
                # kill our current polymer
                c_poly[rxn.pl_idx].pop(cpl_idx)
            # else do nothing
            else:
                react = 0
    # create reaction
    else:
        c_poly[rxn.pl_idx].append(0)
    return react

    # find type of next reaction and fires it
def fire_rxn(rnext, c, t, c_poly, C, T):
    # polymer reaction
    react = 1
    if (rnext.rxn_cat > 0):
        # process
        react = pl_cases(rnext, react, c_poly)        
    if (react):
        c += rnext.stoich()
        T.append(t)
        C.append(c.copy())
    return C, T

def gillespie_nrm(tspan, c_0, c_poly, reactions, dep_graph):
    t = tspan[0]
    c = c_0.copy()
    T = [t]
    C = [c.copy()]

    volume = compute_volume(c_0, c_poly)

    # initialize scheduler
    scheduler = PQDict()

    for rxn in reactions:
        tau = -log(rand())/rxn.propensity(c, volume)
        scheduler[rxn] = t + tau

    # first event
    rnext, tnext = scheduler.topitem()
    t = tnext
    C, T = fire_rxn(rnext, c, t, c_poly, C, T)

    while t < tspan[1]:
        # reschedule dependent reactions
        for rxn in dep_graph[rnext]:
            tau = -log(rand())/rxn.propensity(c, volume)
            scheduler[rxn] = t + tau
        
        # fire the next one!
        rnext, tnext = scheduler.topitem()
        t = tnext
        if (rnext.propensity(c, volume) > 0):
            C, T = fire_rxn(rnext, c, t, c_poly, C, T)
        else:
            print(c_poly)
    return array(T), array(C)

# initialize a polymer species
def poly_init(c, c_poly, sp_idx, poly_list, poly_count):    
    c_poly[poly_dict.get(sp_idx)] = poly_list
    c[sp_idx] = poly_count

def main(argv):
    os.chdir(os.path.dirname(argv[0]) or '.')
    os.chdir('..')
    print(os.getcwd())

    config = configparser.ConfigParser()
    config.read('run.ini')
    run_num = config.getint('setup', 'run_num')
    c_dict = ast.literal_eval(config.get('setup', 'init_config'))
    lead_poly = config.getboolean('setup', 'lead_polymers')
    tspan = (0, config.getfloat('setup', 'tmax'))

    reactions = extract_rxns('data/' + argv[1])
    dep_graph = get_dependencies(reactions)

    # create blank configuration
    c = np.zeros(len(species_list), dtype=int)
    c_poly = [[] for p in poly_dict]

    # initialize the polymers for blank configuration
    for sp in species_list:
        sp_idx = insert_sp(sp)
        if 'bot' in sp:
            if lead_poly is True:
                poly_init(c, c_poly, sp_idx, [0], 1)
            else:
                poly_init(c, c_poly, sp_idx, [], 0)

    # which agents to track and record
    track = config.get('setup', 'track').split(',')
    if len(track) > 1:
        t_idxs = [insert_sp(t.strip()) for t in track]
    else:
        t_idxs = insert_sp(track[0])
    
    fmt = config.get('output', 'format')
    output = config.get('output', 'folder')
    if not os.path.exists(output):
        os.makedirs(output)
    o_name = output + '/' + argv[2]

    # run protocol
    if run_num == 1:
        for sp in c_dict:
            sp_idx = insert_sp(sp)
            if 'bot' in sp:
                poly_init(c, c_poly, sp_idx, c_dict.get(sp), len(c_dict.get(sp)))
            else:
                c[sp_idx] = c_dict.get(sp)
        
        # blanks, remove later
        n = c[insert_sp('I_X')]
        c[insert_sp('B')] = n*n*3 + 40

        startTime = datetime.now()
        T, C = gillespie_nrm(tspan, c, c_poly, reactions, dep_graph)
        print('\nTime elasped: ', datetime.now() - startTime)
        counts = C[:, t_idxs]
        time = T

    else:
        counts = []
        time = []
        n_i = 0
        ci = c.copy() # back-up initial config
        ci_poly = copy.deepcopy(c_poly) # back-up initial polymers
        while n_i < run_num:
            for sp in c_dict:
                sp_idx = insert_sp(sp)
                if isinstance(c_dict.get(sp), list):
                    if 'bot' in sp:
                        poly_init(ci, ci_poly, sp_idx, c_dict.get(sp)[n_i], len(c_dict.get(sp)[n_i]))
                    else:
                        ci[sp_idx] = c_dict.get(sp)[n_i]

                    # blanks, remove later
                    if sp == 'I_X':
                        n = ci[insert_sp('I_X')]
                        ci[insert_sp('B')] = n*n*3 + 40

                else:
                    if 'bot' in sp:
                        poly_init(ci, ci_poly, sp_idx, c_dict.get(sp), len(c_dict.get(sp)))
                    else:
                        ci[sp_idx] = c_dict.get(sp)

            startTime = datetime.now()
            T, C = gillespie_nrm(tspan, ci, ci_poly, reactions, dep_graph)
            print('\nTime elasped: ', datetime.now() - startTime)

            counts.append(C[-1, t_idxs])
            time.append(T[-1])

            ci = c.copy() #restore configuration
            ci_poly = copy.deepcopy(c_poly)

            n_i = n_i + 1
    
    if fmt == 'np' or 'numpy':
        np.save(o_name + '_counts.npy', counts)
        np.save(o_name + '_time.npy', time)        

if __name__ == '__main__':
    main(sys.argv)
