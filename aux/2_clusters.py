from pymol import cmd

cmd.set("sphere_scale", "0.9")

cmd.load("/home/pbarletta/labo/22/migue/run/mut_pdbs/1MHP-H:T33I_H:S52T_H:G53Q_H:G54F.pdb")
cmd.color("salmon", "chain A or chain A")
cmd.color("atomic", "(not elem C)")

cmd.select("id 130+131+136+139+140+141+142+143+144+145+151+2265+2266+2268+2269+2270+2281+2282+2284+2285+2286+2043+2046+2047+2244+742+743+744+745+4685+2318+2319+2320+2321+4690+2322+4692+2323+4691+2757+108+109+111+112+2259+2262+2079+2080+2753+122+2752+705+2300+706+2283+2298+731+748+4664")
cmd.set_name("sele", "cluster_1")
cmd.show("spheres", "cluster_1")
cmd.color("0x9b518f", "cluster_1")
cmd.select("id 4054+1118+4642+4644+1125+4645+4635+4636+4637+4638+4639+1153+1155+4052+4053")
cmd.set_name("sele", "cluster_2")
cmd.show("spheres", "cluster_2")
cmd.color("0x8cddbd", "cluster_2")
