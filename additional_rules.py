# rule to make index file
rule make_index:
    input:
        config["reference"]
    output:
        config["contigNames"]
    script:
        scripts/make_index.py