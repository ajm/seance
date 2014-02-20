parser Newick :
    ignore: " \n\t"
    
    token ID: '[-_\.a-zA-Z0-9]+'
    token NUM: '\-?[0-9]+\.[0-9]+(e\-?[0-9]+)?'

    rule newick: tree ";"   {{ return tree }}
    rule tree: ID ":" NUM   {{ return (ID, float(NUM)) }} 
                | "\\("     {{ result = [] }}
                tree ","    {{ result.append(tree) }}
                tree "\\)"  {{ result.append(tree) }}
                ":" NUM     {{ return (result, float(NUM)) }}
