from plotter import *

path = "Culex_oct23.merged.sorted.locusTag.fixed.gff3"

file = open(path)

charts = dict()
gr = 0
ir = 0

for line in file:
    # Split columns
    cols = line.split("\t")
    # less than 9 columns indicates a header line
    if len(cols) < 9:
        continue
    settings = cols[8].split(";")
    # if it isn't a gene or pseudogene we don't want to graph it
    if cols[2] != "gene" and cols[2] != "pseudogene":
        continue
    # make a list of all settings items that start with "Name="
    name = list(filter(lambda a: a.startswith("Name="), settings))
    # throw an error if there are multiple settings called "Name"
    if len(name) > 1:
        raise ValueError("Expected no more than one name property per line")
    # skip the line if there is no name
    if len(name) == 0: continue
    # [name] -> name
    name = name[0]
    # Remove "Name="
    name = name[5:]
    # Determine color based on whether it's a gr or ir
    if name.lower().startswith("gr"):
        gr += 1
        color = "#b326ff"
    elif name.lower().startswith("ir"):
        ir += 1
        color = "orange"
    else:
        # We don't want to graphs genes other than irs and grs, so skip
        continue
    # Which chromosome is this gene on?
    chromosome = cols[0]
    # Position in Mb
    position = int(cols[3])/1000000
    # Set direction variable to LEFT or RIGHT based on the strandedness
    direction = cols[6]
    if direction == "+":
        direction = RIGHT
    elif direction == "-":
        direction = LEFT
    else:
        raise ValueError("Expected strandedness of \"+\" or \"-\"")
    # Makes a Triangle object and adds it to a list for this chromosome
    if chromosome in charts:
        charts[chromosome].append(Triangle(position, cols[2] == "gene", direction, color))
    else:
        charts[chromosome] = [Triangle(position, cols[2] == "gene", direction, color)]

# Total length for each chromosome
total = {
    "CM027410.1": 132876167,
    "CM027411.1": 225161840,
    "CM027412.1": 201550677
}
# Centromere position for each chromosome
first = {
    "CM027410.1": 63000000,
    "CM027411.1": 109000000,
    "CM027412.1": 95000000
}

# Create output object
output = Output()
# Add graphs in reverse order (3, 2, 1) so when added bottom up it will appear as (1, 2, 3)
charts = sorted(list(charts.items()), key=lambda a: a[0])
charts = charts[::-1]
for i, (chromosome, chart) in enumerate(charts):
    chart.sort(key=lambda a: a.position)
    output.add_graph(Graph(chart, total[chromosome]/1000000, first[chromosome]/1000000, f"Chr. {3-i}"))

print(f"{ir} irs and {gr} grs")
output.finalize(keyinfo={"#b326ff": "GRs", "orange": "IRs"})
output.visualize()