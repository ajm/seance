parser Newick :
    ignore: " \n\t"
    
    token ID: '[-_\.a-zA-Z0-9]+'
    token NUM: '\-?[0-9]+\.[0-9]+(e\-?[0-9]+)?'

    rule newick: tree ";"   {{ return tree }}
    rule tree: ID ":" NUM   {{ return (ID, float(NUM)) }} 
                | "\\("     {{ result = [] ; distance = 0.0 }}
                tree ","    {{ result.append(tree) }}
                tree "\\)"  {{ result.append(tree) }}
                (":" NUM)?  {{ return (result, float(NUM if 'NUM' in locals() else 0.0)) }}

