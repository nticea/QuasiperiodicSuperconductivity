

# Open the file in read mode
file = open("example.txt", "r")

# Read the entire contents of the file
contents = read(file, String)

# Close the file
close(file)

# Define the string to search for
search_string = "pattern"

# Use regular expressions to find matches
matches = matchall(r"\b$search_string\b", contents)

# Print the matches
for match in matches
    println(match.match)
end