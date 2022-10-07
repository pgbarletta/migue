from pymol import cmd

cmd.set("sphere_scale", "0.9")

cmd.load("/home/pbarletta/labo/22/migue/run/mut_pdbs/3NPS-A:R48A.pdb")
cmd.color("salmon", "chain A or chain A")
cmd.color("atomic", "(not elem C)")

cmd.select("id 897+898+899+1954+3205+3206+1955+1956+1958+1959+3203+3204+3216+3220+3189+894+895+896+871+2728+873+874+875+876+877+878+2727+2736+2737+2712+2521+2522+3198+3226+363+3223")
cmd.set_name("sele", "cluster_1")
cmd.show("spheres", "cluster_1")
cmd.color("0x12f3d4", "cluster_1")
cmd.select("id 1345+1346+4590+4597+1336+4603+4604+4605+4607+4608+1333+1334+1335+1338+4606")
cmd.set_name("sele", "cluster_2")
cmd.show("spheres", "cluster_2")
cmd.color("0xcb58a7", "cluster_2")
cmd.select("id 3257+1747")
cmd.set_name("sele", "cluster_3")
cmd.show("spheres", "cluster_3")
cmd.color("0x130dbe", "cluster_3")
cmd.select("id 4855+423")
cmd.set_name("sele", "cluster_4")
cmd.show("spheres", "cluster_4")
cmd.color("0x07c1b5", "cluster_4")
