import sys
import logging
import json
import collections

from copy import deepcopy
from os.path import join

from seance.tools import Pagan
from seance.datatypes import Sequence, IUPAC
from seance.filetypes import FastqFile
from seance.progress import Progress

import dendropy
from dendropy import treecalc


class IdentifiabilityScore(object) :
    def __init__(self) :
        self.log = logging.getLogger('seance')

    def write_fasta(self, seqs, fname) :
        with open(fname, 'w') as f :
            for label,seq in seqs :
                print >> f, ">%s\n%s" % (label, seq)

    def read_nematodes(self, fastq_fname, fprimer, rprimer, diffs, length) :
        tmp = []
        acc2name = {}

        # read in sequences
        fq = FastqFile(fastq_fname)
        fq.open()

        for seq in fq :
            if 'Nematoda' not in seq.id :
                continue
    
            seq.ungap()
            seq.back_translate()

            new_id = seq.id.split()[0][1:]
            tmp.append((new_id, seq.sequence))
    
            acc2name[new_id] = seq.id[seq.id.find('Nematoda'):]

        fq.close()


        # test sequences
        p = Progress("Looking for primer sequences", len(tmp), False)
        p.start()

        tmp2 = []

        for label,seq in tmp :
            findex = IUPAC.seq_position(fprimer, seq, diffs)
 
            if findex != -1 :
                #if IUPAC.seq_position_reverse(rprimer, seq, diffs) != -1 :

                shortseq = seq[findex + len(fprimer) : findex + len(fprimer) + length]
                if 'N' not in shortseq :
                    tmp2.append((label, shortseq))          

            p.increment()

        p.end()

        return tmp2,acc2name

    # XXX   this code is duplicated so many places with tiny tweaks, 
    #       maybe abstract out...
    def distance(self, aligned) :
        leng = float(min(len(aligned[0]), len(aligned[1])))

        last_gap = True
        diff = 0

        for c1,c2 in zip(aligned[0], aligned[1]) :
            if (c1 == '-') and (c2 == '-') :
                continue

            gap = '-' in (c1,c2)

            if last_gap and gap :
                continue

            if not IUPAC.equal(c1, c2) :
                diff += 1

            last_gap = gap

        if last_gap :
            diff -= 1

        return (leng - diff) / leng

    def build_distance_matrix(self, fname) :
        tmp = []
        fq = FastqFile(fname)
        fq.open()

        for seq in fq :
            tmp.append((seq.id[1:], seq.sequence))

        fq.close()


        p = Progress("Calculating distance matrix", (len(tmp)*(len(tmp)-1)) / 2.0, False)
        p.start()

        dist = {}
        for index, labelseq in enumerate(tmp) :
            label1,seq1 = labelseq
            dist[(label1,label1)] = 1.0
            for label2,seq2 in tmp[index+1:] :
                dist[(label1,label2)] = \
                        dist[(label2,label1)] = \
                        self.distance([seq1,seq2])
                
                p.increment()

        p.end()

        return dist, [ label for label, seq in tmp ]

    def get_phylogenetic_metric(self, labels, tree) :
        clustertree = deepcopy(tree)
        clustertree.retain_taxa_with_labels(labels)

        subtree = deepcopy(tree)
        subtree.retain_taxa_with_labels([i.get_node_str() for i in subtree.mrca(taxon_labels=labels).leaf_nodes()])

        return clustertree.length() / subtree.length()

    def assign_scores(self, dist, keys, threshold, tree) :
        scores = []

        p = Progress("Calculating scores", len(keys), False)
        p.start()        

        for k1 in keys :
            tmp = []
            for k2 in keys :
                if dist[(k1,k2)] >= threshold :
                    tmp.append(k2)
            p.increment()
            
            #pdscore = self.get_phylogenetic_metric(tmp, tree)

            scores.append((k1, 1.0, tmp))
            #scores.append((k1, pdscore, tmp))
            #scores.append((k1, 1 / float(len(tmp)), tmp))
        
        p.end()        

        return scores

    def strip_to_genus_nematode(self, s) :
        s = s.split(';')[2:]
        return ';'.join([ i for i in s if i not in ('Nematoda','Enoplea','Chromadorea') ])

    def element_of_list_in_list(self, lista, listb) :
        for i in lista :
            if i in listb :
                return True
        return False

    def merge_taxonomy(self, names) :
        tmp = collections.defaultdict(list)
        for name in names :
            for i,v in enumerate(name.split(';')) :
                tmp[i].append(v)

        s = ""
        for i in sorted(tmp.keys()) :
            if len(set(tmp[i])) == 1 :
                s += (";%s" % tmp[i][0])
            else :
                break

        return s[1:]

    def get_subtree_in_list(self, node, include_list) :
        tmp = []
        for n in node.child_nodes() :
            tmp += self.get_subtree_in_list(n, include_list) 

        name = node.get_node_str()

        if name in include_list :
            tmp.append(name)

        return tmp

    def score(self, silva_fasta, silva_tree, outdir, fprimer, rprimer=None, diffs=2, length=250, threshold=0.99) :
        # read nematodes that the primer would hit, trim to length
        # seqs is list of tuples [(label,seq), (label,seq), ...]
        self.log.info("extracting sequences that match '%s' from '%s' ..." % (fprimer, silva_fasta))
        seqs,acc2name = self.read_nematodes(silva_fasta, fprimer, rprimer, diffs, length)

        # write out nematodes and align using the silva tree
        self.log.info("aligning with pagan...")
        fasta_fname = join(outdir, 'score.fasta')
        self.write_fasta(seqs, fasta_fname)
        align,tree,xml = Pagan().phylogenetic_alignment(fasta_fname, silva_tree)

        # read msa and calculate distance matrix
        self.log.info("building distance matrix and scoring")
        dist,keys = self.build_distance_matrix(align)

        # read tree
        t = dendropy.Tree.get_from_path(silva_tree, schema='newick', as_rooted=True)
        scores = self.assign_scores(dist, keys, threshold, t)

        # XXX THIS IS JUST TRYING TO REDUCE THE DATA FOR TESTING
        include_list = []
        with open('out_denoise_all_040314/newb2.phylogeny.fasta') as f :
            for line in f :
                if line.startswith('>') :
                    tmp = line[1:].split()[0] 
                    if tmp in acc2name :
                        include_list.append(tmp)
        
        tmp = []
        for name,pdscore,cluster in scores :
            if self.element_of_list_in_list(cluster, include_list) :
                tmp.append(name)

        include_list = include_list + tmp

        #include_list = [ i[0] for i in scores ] # ALL!
        # XXX END

        with open('nematode_distances.txt', 'w') as f :
            for k in dist :
                a,b = k
                print >> f, a, b, dist[k]

#        with open('nematode_scores.txt', 'w') as f :
#            print >> f, "score accession name cluster"
#            for name,iscore,cluster in sorted(scores, key=lambda x : x[1], reverse=True) :
#                print >> f, iscore, name, acc2name[name], "(%s)" % ','.join(cluster)

        with open('nematode_scores.txt', 'w') as f :
            print >> f, "accession name taxonomic_cluster_label degree_of_monophyly cluster_size cluster"
            for name,pdscore,cluster in sorted(scores, key=lambda x : x[1], reverse=True) :
                taxonomic_label = self.merge_taxonomy([ acc2name[i] for i in cluster ])
                print >> f, "\t".join([name, acc2name[name], taxonomic_label, str(pdscore), str(len(cluster)), "(%s)" % ','.join(cluster)])


        tmp = []
        for name,pdscore,cluster in scores :
            if name not in include_list :
                continue

            cluster2 = [ i for i in cluster if (i != name) and (i in include_list) ]
            #notcluster = [ i for i in [ j.get_node_str() for j in t.mrca(taxon_labels=cluster2 + [name]).leaf_nodes() ] if (i != name) and (i in include_list) and (i not in cluster2) ]
            notcluster = [ i for i in self.get_subtree_in_list(t.mrca(taxon_labels=cluster2 + [name]), include_list) if (i != name) and (i not in cluster2) ]

            #cluster2 = [ i for i in [ j.get_node_str() for j in t.mrca(taxon_labels=cluster2 + [name]).leaf_nodes() ] if (i != name) and (i in include_list) ]
            #notcluster = []

            tmp.append({
                "name"       : name,
                "label"      : self.strip_to_genus_nematode(acc2name[name]),
                "size"       : pdscore,
                "cluster"    : cluster2,
                "notcluster" : notcluster
                })

        with open('nematode_vis.json', 'w') as f :
            print >> f, json.dumps(tmp)

        tmp = []
        done = {}
        for name,pdscore,cluster in scores :
            if name in include_list :
                for name2 in [ i for i in cluster if (i != name) and (i in include_list) ] :
                    if (name,name2) in done :
                        continue

                    tmp.append({
                        "source" : name,
                        "target" : name2
                        })

                    done[(name,name2)] = done[(name2,name)] = 1

        with open('nematode_links.json', 'w') as f :
            print >> f, json.dumps(tmp)

        # XXX tree things
        t = dendropy.Tree.get_from_path(silva_tree, schema='newick', as_rooted=True)
        t.retain_taxa_with_labels(include_list)
        t.write_to_path('nematode_vis.nwk', 'newick')
        t.ladderize(ascending=False)
        t.write_to_path('nematode_vis_ladder.nwk', 'newick')

        # sequence similarity vs. phylogenetic distance
#        t = dendropy.Tree.get_from_path(silva_tree, schema='newick', as_rooted=True)
#        pdm = treecalc.PatristicDistanceMatrix(t)
#
#        with open('nematode_ss_v_pd.txt', 'w') as f :
#            for i, t1 in enumerate(t.taxon_set):
#                for t2 in t.taxon_set[i+1:] :
#                    if (t1.label,t2.label) not in dist :
#                        continue
#
#                    print >> f, t1.label, t2.label, pdm(t1,t2), dist[(t1.label,t2.label)]

if __name__ == '__main__' :
    IdentifiabilityScore().score('SSU_reference_l5_Silva.fasta', 'SSU_reference_l5_Silva.tree', '.', 'AGRGGTGAAATYCGTGGAC', 'TCTCGCTCGTTATCGGAAT', 2, 250, 0.99)

