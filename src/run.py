from pqdict import PQDict
import numpy as np
from numpy import array, zeros, log, seterr
from numpy.random import rand, randint
from scipy.optimize import curve_fit
from collections import Counter
from datetime import datetime
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

# generation code
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
    pl_idx = poly_dict.get(sp_idx, len(c_poly))
    if (pl_idx == len(c_poly)):
        c_poly.append([])
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
            pl_idx = insert_pl(sp_idx) # as of now only one polymer possible
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

# runtime code
def get_dependencies(reactions):
    dependencies = {}
    for rxn in reactions:
        # dependency list
        dep_list = []
        agent_idx = set(i for i in rxn.count if rxn.count.get(i) != 0)
        # get indices of the reactants/products
        agents = set(species_list[i] for i in agent_idx)
        for rxn2 in reactions:
            agents2 = set(species_list[i] for i in rxn2.r_idx)
            #if bool(agent_idx & set(rxn2.r_idx)):
            #l_agent = i for i in agents2 if i.startswith('L_') and not i.startswith('L_t')
            #if (l_agent == [] and bool(agents & agents2)) or (l_agent != [] and l_agent[0] in agents):
            if bool(agents & agents2) and agents2 & agents != set('B'):
                dep_list.append(rxn2)
        # add key-value pair (reaction and its dependent reactions) to dictionary
        dependencies[rxn] = dep_list
    return dependencies

def compute_volume(c_0, c_poly):
    pl_total = 0
    for pl in c_poly:
        pl_total = pl_total + sum(pl)
    return sum(c_0) + pl_total

def pl_cases(rxn, c_poly, react): #processes effects of polymer lists of polymer reactions
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
                # kill our current stack
                c_poly[rxn.pl_idx].pop(cpl_idx)
            # else do nothing
            else:
                react = 0
    # create reaction
    else:
        c_poly[rxn.pl_idx].append(0)
    return react

    # figures out what type of reaction is next and fires it
def fire_rxn(rnext, c, c_poly, t, C, T):
    # polymer reaction
    react = 1
    if (rnext.rxn_cat > 0):
        # process
        react = pl_cases(rnext, c_poly, react)        
    if (react):
        c += rnext.stoich()
        T.append(t)
        C.append(c.copy())
        #print(species_list[rnext.r_idx[0]] + " + " + species_list[rnext.r_idx[1]] + '->' + \
        #species_list[rnext.p_idx[0]] + " + " + species_list[rnext.p_idx[1]])
        #print(c_poly)
        #print(c[species_list.index('A_X')], c[species_list.index('I_X')])
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
    C, T = fire_rxn(rnext, c, c_poly, t, C, T)

    while t < tspan[1]:
        # reschedule dependent reactions
        for rxn in dep_graph[rnext]:
            tau = -log(rand())/rxn.propensity(c, volume)
            scheduler[rxn] = t + tau
        
        # fire the next one!
        rnext, tnext = scheduler.topitem()
        t = tnext
        if (rnext.propensity(c, volume) > 0):
            C, T = fire_rxn(rnext, c, c_poly, t, C, T)
        else:
            print(c_poly)
    return array(T), array(C)

species_list = []
poly_dict = {}
c_poly = []
reactions = extract_rxns('example_rxns.txt')
dep_graph = get_dependencies(reactions)

n = 100
c_dict = {
    'bot_X': [n],
    'bot_X\'': [0],
    'bot_X\'\'': [0],
    'bot_Y_int': [0],    
    'bot_Y\'': [0],
    'bot_copy': [0],    
    #'X': n,
    'I_X': n,
    'B': n*n*3 + 40,
    'L_1': 1
}
c = np.zeros(len(species_list), dtype=int)
tspan = (0, 1.4e10)

I_Y_idx = insert_sp('I_Y_int')
Y_idx = insert_sp('Y')

for sp in c_dict:
    sp_idx = insert_sp(sp)
    if 'bot' in sp:
        c_poly[poly_dict.get(sp_idx)].extend(c_dict.get(sp))
        c[sp_idx] = len(c_dict.get(sp))
    else:
        c[sp_idx] = c_dict.get(sp)

startTime = datetime.now()

T, C = gillespie_nrm(tspan, c, c_poly, reactions, dep_graph)

print('\nTime elasped: ', datetime.now() - startTime)

#np.save('rseq-pl-counts-100.npy', C[:,I_Y_idx])
#np.save('rseq-sol-counts-100.npy', C[:,Y_idx])
#np.save('rseq-intrt-100.npy', T)

#C0 = np.load('seq-pl-counts-100.npy')
#C1 = np.load('rseq-sol-counts-100.npy')
#T0 = np.load('seq-intrt-100.npy')

C1 = np.load('counts-squaret-100.npy')
T1 = np.load('intrt-squaret-100.npy')

#C1 = np.load('counts-threaded-128.npy')
#T1 = np.load('intrt-threaded-128.npy')

#C1 = np.load('seq-counts.npy')
#T1 = np.load('seq-intrt.npy')


plt.plot(T1, C1, '-b', label = 'threaded')
plt.plot(T, C[:,I_Y_idx], '-r', label = 'sequential') 
plt.plot([T1[-1], T[-1]], [C1[-1], C[-1,I_Y_idx]], '-b')

plt.xlabel('number of interactions')
plt.ylabel('output count (molecules)')    
plt.legend(frameon = True, loc = 'lower right')

#plt.savefig('rseq-{}' .format(n))
plt.savefig('seq-squaret-100')