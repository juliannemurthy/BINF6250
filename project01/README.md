# Introduction
In this project, we are parsing a Variant Call Format (VCF) file to identify and count rare diseases associated with rare genetic variants, based on their ExAC allele frequency. 

# Pseudocode


```


FUNCTION read_file(file_name):

    Create an empty dictionary called disease_counts
    Open file_name

    FOR each line in the file:
        Pass line to parse_line
        Store the returned value as diseases

        FOR each disease in diseases:
            IF disease is already a key in disease_counts:
                Increase its value by 1
            ELSE:
                Add disease as a key with a value of 1

    RETURN disease_counts


#----------------------------------------
#!/usr/bin/env python
from pprint import pprint


# parse_line fucntion
FUNCTION  def parse_line(line: string) -> list of strings

	IF line starts with "#":
			continue

	columns = split line by tab character #.split("\t")

		info_field = columns[7]              
		info_items = split info_field by ";" #.split(";") 

		af_exac_value = NOT FOUND

	FOR each item in info_items:
			IF item starts with "AF_EXAC=": #startswith("AF_EXAC")
			# Split INFO by semicolons to get key-value fields
				raw_value = part of item after "=" #raw_value= item.split("=")[1]
				IF raw_value contains ",":
					af_exac_value = float of (raw_value split by "," )[0]
				ELSE:
					af_exac_value = float(raw_value)
				(BREAK out of loop ) 

		IF af_exac_value == NOT FOUND:
			RETURN empty list

		IF af_exac_value >= 0.0001:
			RETURN empty list       
			
		disease_list = empty list

		FOR each item in info_items:
			IF item starts with "CLNDN=":
				raw_diseases = part of item after "="
				replace "|" with "," in raw_diseases
				disease_names = split raw_diseases by ","

				FOR each name in disease_names:
					trimmed_name = strip whitespace from name #trimmed_name=name.strip()
					IF trimmed_name is empty:
						CONTINUE to next name
					IF trimmed_name == "not_specified" OR trimmed_name == "not_provided":
						CONTINUE to next name
					ADD trimmed_name to disease_list
				BREAK out of loop  # found CLNDN, no need to keep searching

		RETURN disease_list

# -----------------------------------------------------
# read_file function
FUNCTION def read_file(filename: string) -> dictionary
  results = empty dictionary   #results = {}

    OPEN filename FOR READING as file #with open (filename , "r") as file:
        FOR each line in file :
            diseases = parse_line(line)

            FOR each disease in diseases:
                IF disease already a key in results:
                    results[disease] = results[disease] + 1 # results[disease] += 1
                ELSE:
                    results[disease] = 1
    CLOSE file

    RETURN results



if __name__ == "__main__":
    pprint(read_file("clinvar_20190923_short.vcf"))
    




```

# Successes
- This was one of my first times live coding with other people, so getting comfortable talking through the code together was a learning curve at first. It was helpful to work through things I did not know with the group instead of immediately jumping to Google or an LLM. Talking through the logic together made the process feel more collaborative and helped me understand the reasoning behind our choices.
- Having some previous experience with VCF parsing also helped me recognize parts of the file structure and contribute to the discussion more confidently.

# Struggles
- GitHub took the longest to figure out. Cloning, forking, switching branches, pushing changes, and creating pull requests were all new, so we spent a lot of time making sure we were doing each step correctly.
- We are still getting used to reviewing pull requests separately. It would be easier and faster if we could compare all three or more submissions at the same time instead of checking each one individually.
- Since the workflow is still new to us, we were also careful about not accidentally changing the wrong file or working in the wrong branch.
- A lot of the challenge was less about the Python itself and more about making sure everyone’s work was being combined correctly.

# Personal Reflections
## Group Leader
Group leader's reflection on the project

## Other member
## Aamna: 
The biggest learning point for me was understanding the GitHub workflow and how forks, branches, commits, and pull requests all connect. It took some trial and error to make sure I was working in the correct branch and not overwriting anyone else’s work. I also found that reviewing multiple pull requests separately can be a little slow, especially when several people have made changes to the same file. I feel a lot more comfortable with the process now and understand much better how group coding projects are managed through GitHub.

# Generative AI Appendix
As per the syllabus
