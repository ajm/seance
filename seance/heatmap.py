from sys import exit, stderr, argv
from math import log10, pi

from newick import parse as newick
import json

import cairo

# height, width, font size, tree width
# or just box width and everything derives from that

def parse_newick(newick_file) :
    with open(newick_file) as f :
        s = f.read()
    
    return newick('newick', s)

def parse_biom(biom_file) :
    with open(biom_file) as f :
        s = f.read()

    return json.loads(s)

def preprocess_data(tree_data, biom_data) :
    samples = biom_data['columns']
    species = biom_data['rows']
    counts  = biom_data['data']
    
    #
    # we need to remove species and samples that have zero counts
    #
    seen_species, seen_samples = zip(*[ (x,y) for x,y,c in counts ])
    
    seen_species = sorted(set(seen_species))
    seen_samples = sorted(set(seen_samples))

    #
    # get new lists of non-zero sample names and species names
    #
    new_samples = [ samples[i]['id'] for i in seen_samples ]
    new_species = [ species[i]['id'] for i in seen_species ]

    #
    # zero count species need to be pruned from the tree data
    #
    new_tree = prune_tree(tree_data, new_species)

    #
    # resort non-zero species list into order from tree
    #
    tree_order = []
    get_dfs_order(new_tree, tree_order)

#    print >> stderr, "species in count data:", len(seen_species)
#    print >> stderr, "species in tree data:", len(tree_order)

    tmp = {}
    for name,index in zip(new_species, seen_species) :
        tmp[name] = index

    seen_species = [ tmp[i] for i in tree_order ]
    new_species = tree_order

    #
    # change the count data to reflect the ordering from the tree
    #
    not_found = set()
    new_counts = {} # key = (species_index, sample_index)

    for species_x,sample_y,freq in counts :
        try :
            s_x = seen_species.index(species_x)

        except ValueError :
            not_found.add(species[species_x]['id'])
            continue

        s_y = seen_samples.index(sample_y)

        # i want sames on the x and species on the y
        # so swap them
        new_counts[(s_y, s_x)] = freq

    #
    # print out some info
    #
#    print >> stderr, "%d samples, %d species" % (len(samples), len(species))
#    print >> stderr, "%d non-zero samples, %d non-zero species" % (len(new_samples), len(new_species))

    if len(not_found) :
        print >> stderr, "%d samples from the count data were non-zero, but not found in the tree" % (len(not_found))
        for i in not_found :
            print >> stderr, "\t%s" % i

    return { 
             'tree'    : new_tree,
             'species' : species_name_transform(new_species),
             'samples' : new_samples,
             'counts'  : new_counts
            }

def species_name_transform(species) :
    new_species = []

    for s in species :
        tmp = s.split('_')
        identity = "(%s%%)" % tmp[-1]
        new_species.append(' '.join(tmp[1:-1] + [identity]))

    return new_species

def prune_tree(tree, species_names) :
    subtree,distance = tree

    if isinstance(subtree, str) :
        if subtree in species_names :
            #print >> stderr, "prune_tree: found %s" % subtree
            return subtree,distance

        #print >> stderr, "prune_tree: did not find %s" % subtree
        return None,None

    subtree0,distance0 = prune_tree(subtree[0], species_names)
    subtree1,distance1 = prune_tree(subtree[1], species_names)

    if subtree0 :
        if subtree1 :
            return [(subtree1, distance1), (subtree0, distance0)], distance
        else :
            return subtree0,distance+distance0
    else :
        if subtree1 :
            return subtree1,distance+distance1
        else :
            return None,None

def get_dfs_order(tree, tip_list) :
    subtree,distance = tree

    if isinstance(subtree, str) :
        tip_list.append(subtree)
        return

    get_dfs_order(subtree[0], tip_list)
    get_dfs_order(subtree[1], tip_list)

# get depth
def dfs_dimensions(tree, tipcount) :
    subtree,distance = tree
    
    if isinstance(subtree, str) :
        return tipcount + 1, distance
    
    tipcount,d1 = dfs_dimensions(subtree[0], tipcount)
    tipcount,d2 = dfs_dimensions(subtree[1], tipcount)
    
    return tipcount, distance + max(d1, d2)

tree_extent = 0
y_scalar = x_scalar = 1
margin = 10

def dfs_draw(context, tree, tipcount, height) :
    subtree, distance = tree
    distance *= x_scalar

    if isinstance(subtree, str) :
        y = tipcount * y_scalar
        y += (y_scalar / 2.0)
        y += margin

        # draw tip upwards
        context.set_source_rgba(0.0, 0.0, 0.0, 1.0)
        context.move_to(height, y)
        context.line_to(height + distance, y)
        context.stroke()

        context.set_dash([1, 5])
        context.set_line_width(1)
        context.set_source_rgba(0.0, 0.0, 0.0, 0.3)
        context.move_to(tree_extent, y)
        context.line_to(height + distance, y)        
        context.stroke()

        context.set_dash([])
        context.set_line_width(2)

        return tipcount + 1, height, y
    
    new_height = height + distance
    new_tipcount, left_x, left_y = dfs_draw(context, subtree[0], tipcount, new_height)
    new_tipcount, right_x, right_y = dfs_draw(context, subtree[1], new_tipcount, new_height)

    # cross bar 
    context.set_source_rgba(0.0, 0.0, 0.0, 1.0)
    context.move_to(left_x, left_y)
    context.line_to(right_x, right_y)

    # draw line to root
    mid_y = (left_y + right_y) / 2.0
    context.move_to(left_x, mid_y)
    context.line_to(left_x - distance, mid_y)
    
    context.stroke()

    return new_tipcount, left_x - distance, mid_y

def max_text_length(s, context) :
    return max([ context.text_extents(i)[2] for i in s ])

def draw_heatmap(context, count_data, x_start, y_start, block_size) :
    max_count = float(max(count_data.values()))
    max_count = log10(max_count)

    for x,y in count_data :
        context.set_source_rgba(0.0, 0.4, 0.8, log10(count_data[(x,y)]) / max_count)
        
        rect_x = x_start + (x * block_size)
        rect_y = y_start + (y * block_size)

        context.rectangle(rect_x, rect_y, block_size, block_size)
        context.fill()

def draw_species_labels(context, species, x_start, y_start, block_size) :
    context.set_source_rgba(0.0, 0.0, 0.0, 1.0)    

    for ind,s in enumerate(species) :
        tmp = block_size - context.text_extents(s)[3]
        context.move_to(x_start, y_start + ((ind + 1) * block_size) - tmp)
        context.text_path(s)
        context.fill()

def draw_sample_labels(context, samples, x_start, y_start, block_size) :
    context.set_source_rgba(0.0, 0.0, 0.0, 1.0)

    for ind,s in enumerate(samples) :
        tmp = block_size - context.text_extents(s)[3]
        x = x_start + (ind * block_size) + tmp
        y = y_start

        context.translate(x, y)
        context.rotate(pi / 2.0)
        context.translate(-x, -y)

        context.move_to(x, y)
        context.text_path(s)

        context.fill()
        context.identity_matrix()

def heatmap(biomfile, treefile, pdffile, draw_guidelines=False) :
    global x_scalar, y_scalar, margin, tree_extent

    newick_data = parse_newick(treefile)
    biom_data = parse_biom(biomfile)

    data = preprocess_data(newick_data, biom_data)

    # setup cairo with dummy dimensions
    surface = cairo.PDFSurface(pdffile, 0, 0)
    context = cairo.Context(surface)
    
    context.select_font_face("monospace")
    context.set_font_size(8)

    #
    #context.scale(1, 1)
    #context.set_line_width(1)
    #context.set_source_rgb(0, 0, 0)
    #context.set_line_cap(cairo.LINE_CAP_BUTT)
    #context.set_line_join(cairo.LINE_JOIN_ROUND)
    #context.set_dash([10, 5])

    # calculate new variables
    block_len = 10
    tree_blocks = 30
    spacer = 3
    margin = 10

    # calculate based on data
    max_species_label = max_text_length(data['species'], context)
    max_sample_label = max_text_length(data['samples'], context)

    num_species, phylogenetic_height = dfs_dimensions(data['tree'], 0)
    
    heatmap_width = len(data['samples']) * block_len
    heatmap_height = num_species * block_len

    tree_width = tree_blocks * block_len
    tree_height = heatmap_height
    tree_extent = margin + tree_width 

    width =  margin + tree_width + spacer + heatmap_width + spacer + max_species_label + margin
    height = margin + heatmap_height + spacer + max_sample_label + margin

#    print >> stderr, "%d tip nodes in tree" % num_species
#    print >> stderr, "tree height = %f" % height

    x_scalar = tree_width / float(phylogenetic_height)
    y_scalar = tree_height / float(num_species)

#    print >> stderr, "x_scalar = %f" % x_scalar
#    print >> stderr, "y_scalar = %f" % y_scalar


    # change dimensions
    surface.set_size(width, height)

    if draw_guidelines :
        context.set_line_width(1)
        context.set_line_cap(cairo.LINE_CAP_BUTT)
        context.set_line_join(cairo.LINE_JOIN_MITER)
        context.set_source_rgba(1.0, 0, 0, 0.1)
        
        # draw margin
        context.rectangle(margin, margin, width - (margin * 2), height - (margin * 2))

        # draw horizontal guides
        for i in range(1, int(num_species) + 1) :
            context.move_to(margin, margin + (i * block_len)) 
            context.line_to(width - margin, margin + (i * block_len))

        # draw vertical guides
        tmp = margin + tree_width + spacer
        for i in range(len(data['samples']) + 1) :
            context.move_to(tmp + i * block_len, margin)
            context.line_to(tmp + i * block_len, height - margin)

        # spacers
        context.move_to(margin + tree_width, margin)
        context.line_to(margin + tree_width, margin + heatmap_height)
    
        context.move_to(margin + tree_width + spacer + heatmap_width + spacer, margin)
        context.line_to(margin + tree_width + spacer + heatmap_width + spacer, margin + heatmap_height)
    
        context.move_to(margin + tree_width + spacer, margin + heatmap_height + spacer)
        context.line_to(margin + tree_width + spacer + heatmap_width, margin + heatmap_height + spacer)

        context.stroke()

    # draw tree
    context.set_line_width(2)
    context.set_line_cap(cairo.LINE_CAP_BUTT)
    context.set_line_join(cairo.LINE_JOIN_ROUND)    
    context.set_source_rgb(0, 0, 0)

    dfs_draw(context, data['tree'], 0, margin)

    # draw heatmap
    context.set_line_width(1)
    context.set_line_cap(cairo.LINE_CAP_BUTT)
    context.set_line_join(cairo.LINE_JOIN_MITER)

    draw_heatmap(context, data['counts'], margin + tree_width + spacer, margin, block_len)

    # draw gridlines
    # horizontal 
    tmp = margin + tree_width + spacer
    context.set_source_rgba(0.5, 0.5, 0.5, 0.5)
    for i in range(1, int(num_species)) :
        y = margin + (i * block_len)
        context.move_to(tmp, y)
        context.line_to(tmp + heatmap_width, y)

    # vertical 
    for i in range(1, len(data['samples'])) :
        context.move_to(tmp + (i * block_len), margin)
        context.line_to(tmp + (i * block_len), margin + heatmap_height)

    context.rectangle(tmp, margin, heatmap_width, heatmap_height)
    context.stroke()

    draw_species_labels(context, data['species'], margin + tree_width + spacer + heatmap_width + spacer, margin, block_len)

    draw_sample_labels(context, data['samples'], margin + tree_width + spacer, margin + heatmap_height + spacer, block_len)

    context.save()
    surface.finish()

    return 0

def main() :
    if len(argv) != 3 :
        print >> stderr, "Usage: %s <biom file> <tree file>" % argv[0]
        return 1

    biomfile = argv[1]
    treefile = argv[2]

    return heatmap(biomfile, treefile, "tree.pdf")

if __name__ == '__main__' :
    try :
        exit(main())
    except KeyboardInterrupt :
        print >> stderr, "Killed by user"
        exit(1)
