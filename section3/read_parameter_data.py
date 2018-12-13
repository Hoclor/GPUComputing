
# Read the data file
with open('hamilton_parameter_test.dat', 'r') as f:
    # list to store each line entry
    lines = []
    # Read lines until '' is read
    for line in f:
        lines.append(line)

# Split all lines by ',', then remove whitespaces, then convert to floats
processed_lines = []
for index, line in enumerate(lines):
    temp_line = line.split(',')
    for item_index, item in enumerate(temp_line):
        temp_line[item_index] = item.strip(' \n')
        # Convert the first 5 items to ints, the last 2 to floats
        if item_index < 5:
            temp_line[item_index] = int(temp_line[item_index])
        else:
            temp_line[item_index] = float(temp_line[item_index])
    # Remove all lines where line[0] != 1024 as the data for 1536 is incomplete
    if temp_line[0] == 1024:
        processed_lines.append(temp_line[:])
    # Replace the line with the processed line
    lines[index] = temp_line

# Sort the list of lines by the 6th (index=5) item    
processed_lines.sort(key=lambda line: line[5])

# Print the sorted list in reverse
processed_lines.reverse()

for i in processed_lines:
    print(i)

# Un-reverse the list
processed_lines.reverse()