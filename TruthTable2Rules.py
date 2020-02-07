
from pysb.core import MonomerPattern, ComplexPattern, RuleExpression, ReactionPattern, ANY, WILD
from pysb.builder import Builder
from pysb.export.pysb_flat import PysbFlatExporter
from copy import deepcopy
from itertools import combinations, product
import re
from collections import defaultdict


class B_Node:

    def __init__(self):
        self.incidentNodes = []
        self.labels = []
        self.reactions = []
        self.optional_rxns = []
        self.initial = None
        self.motifs = []
        self.boolean = None
        self.boolean_list = None
        self.table = None

class Rxn:

    def __init__(self, mole, react, direct, reacts, targs, rxnTemps):
        self.molecule = mole
        self.reaction = react
        self.direction = direct
        self.reactants = reacts
        self.targets = targs
        self.rxnTemplates = rxnTemps
        self.rxnsParsed = []

def runTT2R(model_name):

    # read in the Booleannet model
    model_nodes = importModel(model_name)

    # read in the logical reaction library
    library = importLibrary('Logical_ODE.mrl')

    MotifBuilder(model_nodes, library)
    # Combine_and_Build(model_nodes, library, input_type='Boolean', output_type='logical', top=1)
    ModelBuilder(0, model_nodes, library, model_name)


def importModel(model_name):

    nodes = {}

    # read in the Boolean model in Booleannet format
    model = open(model_name, 'r')
    bnodes = []
    brules = []
    brules_parced = []
    binputs = []

    within_model = False
    for line in model:

        if '\"\"\"' in line and within_model:
            break

        if within_model:

            if '=' in line and '*' not in line:
                s = line.find(' ')
                nodes[line[:s]] = B_Node()
                if line[line.rfind(' ')+1:-1] == 'False':
                    nodes[line[:s]].initial = 0
                if line[line.rfind(' ')+1:-1] == 'True':
                    nodes[line[:s]].initial = 1

            if '*' in line and ':' in line:
                bnodes.append(line[line.index(':')+2:line.index('*')])
                brules.append('('+line[line.index('=') + 2:-1]+')')
                bp = deepcopy(brules[-1])
                bp = bp.replace('(', ' ( ')
                bp = bp.replace(')', ' ) ')
                bp = bp.strip()
                bp = re.split(r'\s*', bp)
                brules_parced.append(bp)
                binputs.append([])

            if '*' in line and ':' not in line:
                bnodes.append(line[:line.index('*')])
                brules.append(line[line.index('=')+2:-1])
                bp = deepcopy(brules[-1])
                bp = bp.replace('(', ' ( ')
                bp = bp.replace(')', ' ) ')
                bp = bp.strip()
                bp = re.split(r'\s*', bp)
                brules_parced.append(bp)
                binputs.append([])

        if '\"\"\"' in line and not within_model:
            within_model = True

    btable = []

    for i, each in enumerate(brules_parced):
        temp = []
        for item in bnodes:
            if item in each:
                temp.append([each.index(item), item])
        temp.sort()
        for item in temp:
            binputs[i].append(item[1])

    for each in binputs:
        btable.append(list(product([True, False], repeat=len(each))))
    for i, each in enumerate(btable):
        for j, item in enumerate(each):
            btable[i][j] = list(btable[i][j])
    for i, each in enumerate(brules_parced):
        for j, item in enumerate(btable[i]):
            rule = deepcopy(each)
            b_list = deepcopy(item)
            inputs = deepcopy(binputs[i])
            for k, every in enumerate(rule):
                for l, thing in enumerate(inputs):
                    if every == thing:
                        rule[k] = str(b_list[l])
            rule = ' '.join(rule)
            btable[i][j].append(eval(rule))

    for i, each in enumerate(bnodes):
        head = deepcopy(binputs[i])
        head.append(each)
        btable[i].insert(0, head)
    for i, each in enumerate(bnodes):
        nodes[each].boolean = brules[i]
        nodes[each].boolean_list = brules_parced[i]
        nodes[each].table = btable[i]
        for item in binputs[i]:
            nodes[each].incidentNodes.append(item)

    for each in nodes:
        nodes[each].labels.extend(['target', 'reactant'])

    return nodes


def importLibrary(library_name):

    # read in the library of molecules and their associated reactions

    molecule_list = defaultdict(list)

    molecule = None
    reaction = None
    direction = None
    reactants = []
    targets = []
    rxnTemplate = []

    mol_library = open(library_name)
    for line in mol_library:
        if 'molecule:' in line:
            molecule = line.split(':', 1)[1].strip()
        if 'reaction:' in line:
            reaction = line.split(':', 1)[1].strip()
        if 'direction:' in line:
            direction = line.split(':', 1)[1].strip()
        if 'reactant:' in line:
            reactants.append(line.split(':', 1)[1].strip())
        if 'target:' in line:
            targets.append(line.split(':', 1)[1].strip())
        if 'rxnTemplate:' in line:
            rxnTemplate.append(line.split(':', 1)[1].strip())
        if '$$$' in line:
            molecule_list[molecule].append(
                Rxn(molecule, reaction, direction, reactants, targets, rxnTemplate))

            reaction = None
            direction = None
            reactants = []
            targets = []
            rxnTemplate = []

    return molecule_list


class MotifBuilder:

    def __init__(self, nodes, library, hidden=0, max_depth=100):
        self.nodes = nodes
        self.library = library
        self.hidden = hidden
        self.max_depth = max_depth
        self.master = deepcopy(self.nodes)
        self._build_motifs()

    def _findInteractions(self, root_node, target_node, affecting_nodes):

        targ_list = []
        aff_list = []
        interacts = []

        # build lists of reactions for target 'input', affecting 'output', target 'mutual', and affecting 'mutual' nodes
        for target_label in self.master[target_node].labels:
            for mol in self.library[target_label]:
                if mol.direction == 'input':
                    targ_list.append([target_node, mol.reaction, mol.molecule, mol.reactants, 'input'])

        for aff_node in affecting_nodes:
            for aff_label in self.master[aff_node].labels:
                for mol in self.library[aff_label]:
                    if mol.direction == 'output':
                        aff_list.append([aff_node, mol.reaction, mol.molecule, mol.targets, 'output'])

        # group reactions - target 'input', affecting 'output'
        targ_groups = []
        aff_groups = []

        aff_max = 0
        for rxn in aff_list:
            if len(rxn[3]) > aff_max:
                aff_max = len(rxn[3])
        for r in range(aff_max):
            for combos in combinations(targ_list, r+1):
                targ_groups.append(list(combos))

        targ_max = 0
        for rxn in targ_list:
            if len(rxn[3]) > targ_max:
                targ_max = len(rxn[3])
        for r in range(targ_max):
            for combos in combinations(aff_list, r+1):
                aff_groups.append(list(combos))

        # match input and output reactions
        for t_group in targ_groups:
            for a_group in aff_groups:
                match = True

                # match reaction names
                reaction_name = t_group[0][1]
                for t in t_group:
                    if t[1] != reaction_name:
                        match = False
                for a in a_group:
                    if a[1] != reaction_name:
                        match = False

                # match reactants and targets
                if match:
                    for t in t_group:
                        t3 = deepcopy(t[3])
                        if len(t3) != len(a_group):
                            match = False
                        a2 = []
                        for a in a_group:
                            a2.append(a[2])
                        if set(t3) != set(a2):
                            match = False
                    for a in a_group:
                        a3 = deepcopy(a[3])
                        if len(a3) != len(t_group):
                            match = False
                        t2 = []
                        for t in t_group:
                            t2.append(t[2])
                        if set(a3) != set(t2):
                            match = False

                if match:
                    interact = [[], [], [], [], 'd', reaction_name]
                    for t in t_group:
                        interact[1].append(t[0])
                        interact[3].append(t[2])
                    for a in a_group:
                        interact[0].append(a[0])
                        interact[2].append(a[2])
                    interacts.append(interact)

        return interacts

    def _build_motifs(self):

        for node in self.nodes:
            nodelist = deepcopy(self.nodes[node].incidentNodes)
            self.nodes[node].motifs.append([[node], None, None, None, None, None, None, None])
            interactions = []
            interactions += self._findInteractions(node, node, nodelist)
            for each in interactions:
                self.nodes[node].motifs.append(each)


class ModelBuilder(Builder):
    """
    Assemble a PySB model from a Boolean model.
    """
    def __init__(self, num, nodes, library, model_name, iv_type='distribute'):

        super(ModelBuilder, self).__init__()
        self.num = num
        self.nodes = nodes
        self.library = library
        self.model_name = model_name
        self.monomer_info = defaultdict(list)
        self.action_info = defaultdict(list)
        self.base_states = defaultdict(list)
        self.active_states = defaultdict(list)
        self.inactive_states = defaultdict(list)
        self.iv_inactive_states = defaultdict(list)
        self.iv_type = iv_type
        self._build()
        self._export()

    def _build(self):

        self._parse()
        self._get_monomer_info()
        self._add_monomers()
        self._add_rules()
        self._add_initials()
        self._add_observables()

    def _parse(self):

        # parse library reactions
        for every in self.library:
            for thing in self.library[every]:
                for each in thing.rxnTemplates:
                    mols = re.split(r'\s*>>\s*|\s*\+\s*|\s*<>\s*|\s*%\s*', each)
                    ops = re.findall(r'\s*>>\s*|\s*\+\s*|\s*<>\s*|\s*%\s*', each)
                    parced = []
                    for m in mols:
                        parced.append([])
                        sites = []
                        states = []
                        if '(' in m:
                            parced[-1].append(m[:m.index('(')])
                            if '()' not in m:
                                ms = re.split(r'\s*=\s*|\s*,\s*', m[m.index('(') + 1:-1])
                                sites.extend(ms[::2])
                                for s in ms[1::2]:
                                    states.append(s)
                        else:
                            parced[-1].append(m)
                        parced[-1].append(sites)
                        parced[-1].append(states)
                    parced_rxn = [parced.pop(0)]
                    for i, m in enumerate(parced):
                        parced_rxn.append(ops[i])
                        parced_rxn.append(m)
                    thing.rxnsParsed.append(parced_rxn)

        for node in self.nodes:
            for inter in self.nodes[node].motifs:
                if inter[1]:
                    for every in self.library[inter[3][0]]:
                        if every.reaction == inter[5]:
                            pass

    def _export(self):

        f = open(self.model_name.split('.')[0] + '_pysb' + '.py', 'w+')
        f.write(PysbFlatExporter(self.model).export())
        f.close()

        f = open(self.model_name.split('.')[0] + '_pysb' + '.py', "r")
        contents = f.readlines()
        f.close()

        contents.insert(2, 'import numpy as np\nfrom pysb.simulator import ScipyOdeSimulator\nimport pylab as pl\nimport matplotlib.pyplot as plt\n')

        f = open(self.model_name.split('.')[0] + '_pysb' + '.py', "w")
        contents = "".join(contents)
        f.write(contents)
        f.close()


        # f = open('model_' + str(self.num) + '.py', 'w+')
        # f.write(PysbFlatExporter(self.model).export())
        # f.close()
        #
        # f = open('model_' + str(self.num) + '.py', "r")
        # contents = f.readlines()
        # f.close()
        #
        # contents.insert(2, 'import numpy as np\nfrom pysb.integrate import Solver\nimport pylab as pl\nimport matplotlib.pyplot as plt\n')
        #
        # f = open('model_' + str(self.num) + '.py', "w")
        # contents = "".join(contents)
        # f.write(contents)
        # f.close()

    def _find_sites(self, interaction):

        if interaction[1]:

            monomer_names = []
            monomer_labels = []
            monomer_sites = []

            # get names and labels from interaction and initialize lists
            for i, each in enumerate(interaction[0]):
                monomer_names.append(each)
                monomer_labels.append(interaction[2][i])
                monomer_sites.append({})
            for i, each in enumerate(interaction[1]):
                monomer_names.append(each)
                monomer_labels.append(interaction[3][i])
                monomer_sites.append({})

            # break down reaction template and add sites and states
            for every in self.library[interaction[3][0]]:
                if every.reaction == interaction[5]:
                    for thing in every.rxnTemplates:
                        rxn_template = thing
                        rxnTemp = re.split(r'\s*:', rxn_template)[0]
                        mol_list = re.split(r'\s*>>\s*|\s*\+\s*|\s*<>\s*|\s*%\s*', rxnTemp)
                        for m in mol_list:
                            if '(' in m and '()' not in m:
                                mol = m[:m.index('(')]
                                sites = re.split(r'\s*=\s*|\s*,\s*', m[m.index('(') + 1:-1])[::2]
                                states = re.split(r'\s*=\s*|\s*,\s*', m[m.index('(') + 1:-1])[1::2]
                                for i, each in enumerate(sites):
                                    if each not in monomer_sites[monomer_labels.index(mol)]:
                                        monomer_sites[monomer_labels.index(mol)][each] = []
                                        monomer_sites[monomer_labels.index(mol)][each].append(states[i])
                                    else:
                                        if states[i] not in monomer_sites[monomer_labels.index(mol)][each]:
                                            monomer_sites[monomer_labels.index(mol)][each].append(states[i])

            # rename the sites appropriately
            for i, each in enumerate(monomer_names):
                for item in monomer_sites[i]:
                    if item in monomer_labels:
                        monomer_sites[i][monomer_names[monomer_labels.index(item)]] = monomer_sites[i].pop(item)

            # add sites to monomer_info
            for i, each in enumerate(monomer_names):
                for item in monomer_sites[i]:
                    if item not in self.monomer_info[each][0]:
                        self.monomer_info[each][0].append(item)

            # initiate state sites in monomer_info
            for i, each in enumerate(monomer_names):
                for j, item in enumerate(monomer_sites[i]):
                        self.monomer_info[each][1][item] = []

            # add state sites to monomer_info
            for i, each in enumerate(monomer_names):
                for j, item in enumerate(monomer_sites[i]):
                    for every in monomer_sites[i][item]:
                        if every not in self.monomer_info[each][1][item]:
                            self.monomer_info[each][1][item].append(every)

    def _get_monomer_info(self):

        # initiate monomer information
        for motif in self.nodes:
            for interaction in self.nodes[motif].motifs:
                if interaction[1]:
                    for i, species in enumerate(interaction[0]):
                        if species not in self.monomer_info:
                            self.monomer_info[species] = [[], {}]
                    for i, species in enumerate(interaction[1]):
                        if species not in self.monomer_info:
                            self.monomer_info[species] = [[], {}]
                else:
                    for i, species in enumerate(interaction[0]):
                        if species not in self.monomer_info:
                            self.monomer_info[species] = [[], {}]

        # find binding and state sites; inferred from motifs/self interactions and library

        for node in self.nodes:
            for interaction in self.nodes[node].motifs:
                self._find_sites(interaction)

    def _add_monomers(self):

        for each in self.monomer_info:
            self.monomer(each, self.monomer_info[each][0], self.monomer_info[each][1])

    def _get_action_info(self):

        for motif in self.nodes:
            motif_action = []
            for interaction in self.nodes[motif].motifs:
                if interaction[1]:
                    reactants = interaction[0]
                    targets = interaction[1]
                    target_states = [[[], [], []] for _ in range(len(targets))]
                    for every in self.library[interaction[3][0]]:
                        if every.reaction == interaction[5]:
                            for thing in every.rxnTemplates:
                                mol_list = re.split(r'\s*>>\s*|\s*\+\s*|\s*<>\s*|\s*%\s*', thing)
                                mol_list2 = []
                                for m in mol_list:
                                    mol_list2.append(re.split(r'\s*\(\s*|\s*,\s*|\s*=\s*', m[:-1]))
                                for i, each in enumerate(mol_list2):
                                    for j, item in enumerate(each):
                                        for k, stuff in enumerate(interaction[2]):
                                            if stuff == item:
                                                mol_list2[i][j] = interaction[0][k]
                                            if stuff + '_s' == item:
                                                mol_list2[i][j] = interaction[0][k] + '_s'
                                        for k, stuff in enumerate(interaction[3]):
                                            if stuff == item:
                                                mol_list2[i][j] = interaction[1][k]
                                            if stuff + '_s' == item:
                                                mol_list2[i][j] = interaction[1][k] + '_s'
                                reactant_present = [False for _ in range(len(reactants))]
                                for each in mol_list2:
                                    for i, item in enumerate(reactants):
                                        if each[0] == item:
                                            reactant_present[i] = True
                                if all(reactant_present):
                                    for i, each in enumerate(mol_list2):
                                        if len(each) == 2 and each[1] == '':
                                            mol_list2[i].append('')
                                    for each in mol_list2:
                                        for i, item in enumerate(targets):
                                            if each[0] == item and each[1]:
                                                if not target_states[i][0]:
                                                    target_states[i][0] = each[1::2]
                                                    target_states[i][1] = each[2::2]
                                                else:
                                                    target_states[i][2] = each[2::2]
                                    for i in range(len(target_states)):
                                        if not target_states[i][2]:
                                            target_states[i][2] = target_states[i][1]
                    motif_action.append([interaction[0], interaction[1], target_states])
            self.action_info[motif] = motif_action

    def _add_rules(self):

        self._get_action_info()

        # create dictionary of base states
        for each in self.monomer_info:
            mon_obj = self.model.monomers[each]
            states = {}
            for item in mon_obj.sites:
                if item in mon_obj.site_states:
                    states[item] = 'off'
                else:
                    states[item] = 'None'
            state_list1 = []
            state_list2 = []
            for item in states:
                state_list1.append(item)
                state_list2.append(states[item])
            state_list3 = [deepcopy(state_list1), deepcopy(state_list2)]
            self.base_states[each] = state_list3

        # create dictionary of active states based on the Boolean equations
        for each in self.monomer_info:
            a_states = []
            i_states = []

            if each in self.nodes:  # if model node; else proxy node
                table = self.nodes[each].table
                for item in table[1:]:  # find all combinations of active incident nodes that activate the target
                    if item[-1]:
                        boole = defaultdict(bool)  # dictionary of incident node Boolean states on a given row
                        for i, every in enumerate(table[0][:-1]):
                            boole[every] = item[i]

                        # run through all actions for a given motif; consider only those acting on the target node
                        # add sites for active incident nodes
                        sites = {}
                        for every in self.action_info[each]:
                            if every[1][0] == each:
                                for thing in every[0]:
                                    if thing in boole and boole[thing]:
                                        for j, lotsa in enumerate(every[2][0][0]):
                                            if lotsa not in sites and every[2][0][1][j] != every[2][0][2][j]:
                                                sites[lotsa] = every[2][0][2][j]

                        # add sites for inactive incident nodes
                        for every in self.action_info[each]:
                            if every[1][0] == each:
                                for thing in every[0]:
                                    if thing in boole and not boole[thing]:
                                        for j, lotsa in enumerate(every[2][0][0]):
                                            if lotsa not in sites and every[2][0][1][j] != every[2][0][2][j]:
                                                sites[lotsa] = every[2][0][1][j]

                        if sites not in a_states:
                            a_states.append(sites)

                    else:
                        boole = defaultdict(bool)  # dictionary of incident node Boolean states on a given row
                        for i, every in enumerate(table[0][:-1]):
                            boole[every] = item[i]

                        # run through all actions for a given motif; consider only those acting on the target node
                        # add sites for active incident nodes
                        sites = {}
                        for every in self.action_info[each]:
                            if every[1][0] == each:
                                for thing in every[0]:
                                    if thing in boole and boole[thing]:
                                        for j, lotsa in enumerate(every[2][0][0]):
                                            if lotsa not in sites and every[2][0][1][j] != every[2][0][2][j]:
                                                sites[lotsa] = every[2][0][2][j]

                        # add sites for inactive incident nodes
                        for every in self.action_info[each]:
                            if every[1][0] == each:
                                for thing in every[0]:
                                    if thing in boole and not boole[thing]:
                                        for j, lotsa in enumerate(every[2][0][0]):
                                            if lotsa not in sites and every[2][0][1][j] != every[2][0][2][j]:
                                                sites[lotsa] = every[2][0][1][j]
                        if sites not in i_states:
                            i_states.append(sites)

                    self.active_states[each] = a_states
                    self.inactive_states[each] = i_states

        # create rules
        used_interactions = []
        for node in self.nodes:
            for interaction in self.nodes[node].motifs:
                if interaction[1] and interaction not in used_interactions:
                    used_interactions.append(interaction)

                    # retrieve rxn, reaction template, and instructions
                    current_rxn = None
                    for rxn in self.library[interaction[3][0]]:
                        if rxn.reaction == interaction[5]:
                            current_rxn = rxn

                    n = 0
                    for temp in current_rxn.rxnTemplates:

                        rxnTemp1 = re.split(r'\s*:\s*', temp)
                        rxnTemp = rxnTemp1[0]

                        # split template
                        rxn_split = re.split(r'\s*>>\s*|\s*\+\s*|\s*<>\s*|\s*%\s*', rxnTemp)
                        rxn_split_mols = []
                        for each in rxn_split:
                            mol = each.split('(')[0]
                            sites = each.split('(')[1][:-1]
                            site_names = []
                            site_values = []
                            if sites:
                                sites = sites.split(',')
                                for item in sites:
                                    site_names.append(item.split('=')[0].strip())
                                    site_values.append(item.split('=')[1].strip())
                            rxn_split_mols.append([mol, site_names, site_values])
                        ops = re.findall(r'\s*>>\s*|\s*\+\s*|\s*<>\s*|\s*%\s*', rxnTemp)
                        for i, each in enumerate(ops):
                            ops[i] = ops[i].strip()
                        rxn_split_parsed = [rxn_split_mols[0]]
                        for i, each in enumerate(ops):
                            rxn_split_parsed.append(each)
                            rxn_split_parsed.append(rxn_split_mols[i + 1])

                        # Here we make substitutions from the interaction to the template.
                        # Note that the motif generation step does not differentiate between identical
                        # species from different cells. Because the substitutions relies on information
                        # from that step we must account for those identical species here.
                        # !!!!!!!!!!!!!! NEEDS SIMPLIFICATION !!!!!!!!!!!!!!!!!!!

                        for j, each in enumerate(interaction[2]):
                            for i, item in enumerate(rxn_split_parsed):
                                if item != '>>' and item != '<>' and item != '+' and item != '%':
                                    if item[0] == each:
                                        rxn_split_parsed[i][0] = interaction[0][j]
                                    for k, every in enumerate(item[1]):
                                        if every == each:
                                            rxn_split_parsed[i][1][k] = interaction[0][j]
                                        if every[:every.rfind('_')] == each and (every[every.rfind('_')+1:].isdigit() or every[every.rfind('_')+1:] == 's'):
                                            rxn_split_parsed[i][1][k] = interaction[0][j] + every[every.rfind('_'):]

                        for j, each in enumerate(interaction[3]):
                            for i, item in enumerate(rxn_split_parsed):
                                if item != '>>' and item != '<>' and item != '+' and item != '%':
                                    if item[0] == each:
                                        rxn_split_parsed[i][0] = interaction[1][j]
                                    for k, every in enumerate(item[1]):
                                        if every == each:
                                            rxn_split_parsed[i][1][k] = interaction[1][j]
                                        if every[:every.rfind('_')] == each and (every[every.rfind('_')+1:].isdigit() or every[every.rfind('_')+1:] == 's'):
                                            rxn_split_parsed[i][1][k] = interaction[0][j] + every[every.rfind('_'):]

                        rxn_split_parsed_reverse = deepcopy(rxn_split_parsed)
                        rxn_split_parsed_reverse[2][2], rxn_split_parsed_reverse[6][2] = rxn_split_parsed_reverse[6][2], rxn_split_parsed_reverse[2][2]

                        rsp_list = []
                        react = rxn_split_parsed[0][0]
                        for each in self.active_states[react]:
                            rsp = deepcopy(rxn_split_parsed)
                            for item in each:
                                rsp[0][1].append(item)
                                rsp[0][2].append(each[item])
                                rsp[4][1].append(item)
                                rsp[4][2].append(each[item])
                            if rsp not in rsp_list:
                                rsp_list.append(rsp)

                        for each in self.inactive_states[react]:
                            rspr = deepcopy(rxn_split_parsed_reverse)
                            for item in each:
                                rspr[0][1].append(item)
                                rspr[0][2].append(each[item])
                                rspr[4][1].append(item)
                                rspr[4][2].append(each[item])
                            if rspr not in rsp_list:
                                rsp_list.append(rspr)

                        for rsp in rsp_list:

                            # define the rule rule_name
                            rule_name = ''
                            for item in interaction[0]:
                                if item:
                                    rule_name += item + '_'
                                else:
                                    rule_name = rule_name[:-1]
                            rule_name += interaction[5] + '_'
                            for item in interaction[1]:
                                if item:
                                    rule_name += item + '_'
                            rule_name += str(n)

                            # define monomer patterns
                            mon_pats = []
                            for item in rsp:
                                if item == '+' or item == '%' or item == '>>' or item == '<>':
                                    mon_pats.append(item)
                                else:
                                    if item[0] == 'None':
                                        mon_pats.append('None')
                                    else:
                                        mon_states = {}
                                        for i, every in enumerate(item[1]):
                                            mon_states[every] = item[2][i]
                                        mon_obj = self.model.monomers[item[0]]
                                        mon_pats.append(MonomerPattern(mon_obj, mon_states, None))

                            # define complex patterns
                            com_pats_temp = [[]]
                            for item in mon_pats:
                                if item == '>>' or item == '<>':
                                    com_pats_temp.extend([item, []])
                                elif item == '+':
                                    com_pats_temp.append([])
                                elif item == '%':
                                    pass
                                else:
                                    com_pats_temp[-1].append(item)
                            com_pats = []
                            for item in com_pats_temp:
                                if item == '>>' or item == '<>':
                                    com_pats.append(item)
                                elif item == ['None']:
                                    pass
                                else:
                                    com_pats.append(ComplexPattern(item, None))

                            # define reversibility and split patterns into reactants and products
                            react_com_pats = []
                            prod_com_pats = []
                            carrot = 0
                            reversible = None
                            for item in com_pats:
                                if item == '<>':
                                    carrot = 1
                                    reversible = True
                                elif item == '>>':
                                    carrot = 1
                                    reversible = False
                                else:
                                    if carrot == 0:
                                        react_com_pats.append(item)
                                    if carrot == 1:
                                        prod_com_pats.append(item)
                            order = [len(react_com_pats), len(prod_com_pats)]

                            # define rule expression
                            rule_exp = RuleExpression(ReactionPattern(react_com_pats),
                                                      ReactionPattern(prod_com_pats),
                                                      reversible)

                            # create rule
                            if reversible:
                                forward = rule_name + '_' + str(order[0]) + 'kf'
                                self.parameter(forward, 1)
                                reverse = rule_name + '_' + str(order[1]) + 'kr'
                                self.parameter(reverse, 1)
                                self.rule(rule_name, rule_exp, self.model.parameters[forward],
                                          self.model.parameters[reverse])
                            else:
                                forward = rule_name + '_' + str(order[0]) + 'kc'
                                self.parameter(forward, 1)
                                self.rule(rule_name, rule_exp, self.model.parameters[forward])
                            n += 1

    def _add_initials(self):

        if self.iv_type == 'calibrate':

            for each in self.monomer_info:

                mon_obj = self.model.monomers[each]

                state_list = [self.model.monomers[each].sites]
                for item in list(product(['on', 'off'], repeat=len(self.model.monomers[each].sites))):
                    state_list.append(list(item))

                for j, every in enumerate(state_list):

                    if j > 0:
                        index = j-1

                        init_name = each
                        active_state = {}
                        for k, thing in enumerate(mon_obj.sites):
                            active_state[thing] = every[k]
                        init_name += '_' + str(index) + '_0'
                        self.parameter(init_name, 1)
                        self.initial(MonomerPattern(mon_obj, active_state, None), self.model.parameters[init_name])

        if self.iv_type == 'pick':

            # This method picks the state that is most different from all opposing states

            for each in self.monomer_info:
                current_best = [0, None]
                if self.nodes[each].initial:
                    for item in self.active_states[each]:
                        item_score = 100
                        for every in self.inactive_states[each]:
                            every_score = 0
                            for i, thing in enumerate(item):
                                if item[thing] != every[thing]:
                                    every_score += 1
                            if every_score < item_score:
                                item_score = deepcopy(every_score)
                        if item_score > current_best[0]:
                            current_best = [deepcopy(item_score), item]
                else:
                    for item in self.inactive_states[each]:
                        item_score = 100
                        for every in self.active_states[each]:
                            every_score = 0
                            for i, thing in enumerate(item):
                                if item[thing] != every[thing]:
                                    every_score += 1
                            if every_score < item_score:
                                item_score = deepcopy(every_score)
                        if item_score > current_best[0]:
                            current_best = [deepcopy(item_score), item]
                mon_obj = self.model.monomers[each]
                init_name = each
                init_state = deepcopy(current_best[1])
                for every in mon_obj.sites:
                    if every not in init_state:  # absence of site_states is assumed
                        init_state[every] = None
                init_name += '_0'
                self.parameter(init_name, 1)
                self.initial(MonomerPattern(mon_obj, init_state, None), self.model.parameters[init_name])

        if self.iv_type == 'distribute':

            # This method distributes '1' over the active or inactive states

            for each in self.monomer_info:

                mon_obj = self.model.monomers[each]

                state_list = [self.model.monomers[each].sites]
                for item in list(product(['on', 'off'], repeat=len(self.model.monomers[each].sites))):
                    state_list.append(list(item))

                # assume initial boolean values are given
                value = float(self.nodes[each].initial) / len(self.active_states[each])
                for i, item in enumerate(self.active_states[each]):

                    index = None
                    for j, every in enumerate(state_list):
                        if j > 0:
                            match = True
                            for k, thing in enumerate(every):
                                if item[state_list[0][k]] != thing:
                                    match = False
                            if match:
                                index = j-1

                    init_name = each
                    active_state = deepcopy(item)
                    for every in mon_obj.sites:
                        if every not in active_state:
                            active_state[every] = None
                    init_name += '_' + str(index) + '_0'
                    self.parameter(init_name, value)
                    self.initial(MonomerPattern(mon_obj, active_state, None), self.model.parameters[init_name])

                value = (1.0 - float(self.nodes[each].initial)) / len(self.inactive_states[each])
                for i, item in enumerate(self.inactive_states[each]):

                    index = None
                    for j, every in enumerate(state_list):
                        if j > 0:
                            match = True
                            for k, thing in enumerate(every):
                                if item[state_list[0][k]] != thing:
                                    match = False
                            if match:
                                index = j - 1

                    init_name = each
                    inactive_state = deepcopy(item)
                    for every in mon_obj.sites:
                        if every not in inactive_state:
                            inactive_state[every] = None
                    init_name += '_' + str(index) + '_0'
                    self.parameter(init_name, value)
                    self.initial(MonomerPattern(mon_obj, inactive_state, None), self.model.parameters[init_name])

    def _add_observables(self):

        # add all active states combinations to observables
        # for each in self.active_states:
        #     for i, item in enumerate(self.active_states[each]):
        #
        #         obs_name = each + '_' + str(i) + '_obs'
        #
        #         mon_states = {}
        #         for every in item:
        #             for thing in self.action_info[each]:
        #                 for lotsa in thing[2][0][0]:
        #                     if lotsa == every:
        #                         mon_states[every] = item[every]
        #         for every in mon_states:
        #             if mon_states[every].isdigit():
        #                 mon_states[every] = ANY
        #             if mon_states[every] == 'None':
        #                 mon_states[every] = None
        #         mon_pat = MonomerPattern(self.model.monomers[each], mon_states, None)
        #
        #         self.observable(obs_name, mon_pat)

        # adds all possible site combinations to observables
        for each in self.monomer_info:

            state_list = [self.model.monomers[each].sites]
            for item in list(product(['on', 'off'], repeat=len(self.model.monomers[each].sites))):
                state_list.append(list(item))

            for i, every in enumerate(state_list):
                if i > 0:

                    obs_name = each + '_' + str(i-1) + '_obs'
                    mon_states = {}
                    for j, thing in enumerate(every):
                        mon_states[state_list[0][j]] = thing
                    pat = MonomerPattern(self.model.monomers[each], mon_states, None)
                    self.observable(obs_name, pat)

        # adds only monomers to observables
        for each in self.monomer_info:
            obs_name = each + '_obs'
            self.observable(obs_name, self.model.monomers[each])
