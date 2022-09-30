from pymol import cmd

cmd.set("sphere_scale", "0.9")

cmd.load("/home/pbarletta/labo/22/migue/run/mut_pdbs/1BJ1-W:H86A.pdb")
cmd.color("salmon", "chain H or chain L or chain L")
cmd.color("atomic", "(not elem C)")

cmd.select("id 259+260+269+1038+1039+1040+273+274+275+1044+1045+1046+5654+540+5671+5672+5677+942+5678+5680+945+5681+947+948+2936+5682+5689+5690+5695+5696+450+967+5703+969+455+5707+5708+5709+456+975+5712+977+5714+979+976+981+978+5719+980+5721+5722+5723+5724+991+5728+5729+993+2912+5734+5735+5736+5739+5745+2934+1015+457+5628+5629+5632+258+492+493+495+5621+494+496+5331+5332+5333+5698+491+5301+1016")
cmd.set_name("sele", "cluster_1")
cmd.show("spheres", "cluster_1")
cmd.color("0x495f6a", "cluster_1")
