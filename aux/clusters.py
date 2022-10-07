from pymol import cmd

cmd.set("sphere_scale", "0.9")

cmd.load("/home/pbarletta/labo/22/migue/run/mut_pdbs/1BJ1-W:H86A.pdb")
cmd.color("salmon", "chain H or chain L or chain L")
cmd.color("atomic", "(not elem C)")

cmd.select("id 1038+1039+1040+1044+1045+1046+5654+5672+5677+5678+5680+945+5681+947+5682+5689+5690+5695+5696+450+455+456+457+5703+969+5709+975+5712+976+5714+977+979+978+980+5719+981+5724+991+2912+5729+5728+993+5734+5735+5736+2934+2936+5628+5629+5632+5707+5708+269+942+273+275+494+496+5331+5332+5333+492+5621+493+495+259+5723+5721+260+5722+967")
cmd.set_name("sele", "cluster_1")
cmd.show("spheres", "cluster_1")
cmd.color("0x89c980", "cluster_1")
