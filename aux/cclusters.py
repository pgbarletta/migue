from pymol import cmd

cmd.set("sphere_scale", "0.9")

cmd.load("/home/pbarletta/labo/22/migue/data/AB-Bind-Database-master/1DVF.pdb")
cmd.color("salmon", "chain A or chain B or chain B")
cmd.color("atomic", "(not elem C)")

cmd.select("id 803+804+806+784+785+786+787+788+792+382+383+384+385+387+388+390+391+413+412+808+807")
cmd.set_name("sele", "cluster_1")
cmd.show("spheres", "cluster_1")
cmd.color("0xe89542", "cluster_1")
cmd.select("id 710+711+713+714+719+209+823")
cmd.set_name("sele", "cluster_2")
cmd.show("spheres", "cluster_2")
cmd.color("0xf622c4", "cluster_2")
cmd.select("id 209+775")
cmd.set_name("sele", "cluster_3")
cmd.show("spheres", "cluster_3")
cmd.color("0x3b2e83", "cluster_3")
