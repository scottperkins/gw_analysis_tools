import os 

headerList = os.listdir("../../include/gwat/")
for h in headerList:
    if h[-2:] == ".h":
        print("Generating rst file for {}".format(h))
        with open(h[:-2]+".rst","w") as f:
            f.write(".. _api_{}:\n\n".format(h[:-2].lower()))
            f.write("{}\n".format(h[:-2]))
            f.write("="*len(h[:-2])+"\n\n")
            f.write(".. doxygenfile:: {}\n".format(h))
            f.write("\t:project: GW Analysis Tools\n".format(h))
