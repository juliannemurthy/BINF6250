# Introduction
In this project, we are parsing a Variant Call Format (VCF) file to identify and count rare diseases associated with rare genetic variants, based on their ExAC allele frequency. 

# Pseudocode


```
#!/usr/bin/env python
from pprint import pprint


# parse_line fucntion
FUNCTION parse_line(line: string) -> list of strings

	IF line starts with "#":
		RETURN empty list

	columns = split line by tab character #.split("\t")
	info_field = columns[7]              
	info_items = split info_field by ";" #.split(";") 

	af_exac_value = NOT FOUND

	FOR each item in info_items:
		IF item starts with "AF_EXAC=": #startswith("AF_EXAC")
			raw_value = part of item after "="
            af_exac_value = float of raw_value
			(BREAK out of loop ) 

	IF af_exac_value == NOT FOUND:
		RETURN empty list

	IF af_exac_value >= 0.0001:
		RETURN empty list       
		
	disease_list = empty list

	FOR each item in info_items:
		IF item starts with "CLNDN=":
			raw_diseases = part of item after "="
			disease_names = split raw_diseases by "|"
			FOR each name in disease_names:
				trimmed_name = strip whitespace from name #trimmed_name=name.strip()
				IF trimmed_name is empty:
					CONTINUE to next name
				IF trimmed_name == "not_specified" OR trimmed_name == "not_provided":
					CONTINUE to next name
				ADD trimmed_name to disease_list
			BREAK out of loop  

	RETURN disease_list

# -----------------------------------------------------
# read_file function
FUNCTION read_file(filename: string) -> dictionary
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
- We were succesfully able to navigate how to run Git both in the command line and in the browser, gaining confidence in the process.

# Struggles
- GitHub took the longest to figure out. Cloning, forking, switching branches, pushing changes, and creating pull requests were all new, so we spent a lot of time making sure we were doing each step correctly.
- We are still getting used to reviewing pull requests separately. It would be easier and faster if we could compare all three or more submissions at the same time instead of checking each one individually.
- Since the workflow is still new to us, we were also careful about not accidentally changing the wrong file or working in the wrong branch.
- A lot of the challenge was less about the Python itself and more about making sure everyone’s work was being combined correctly.

# Personal Reflections
## Group Leader
### Yulia
As the group leader, it was a valuable learning experience to set up the Github repo with branches, learn about pull requests, and navigate how to accept, reject, or merge edits. My other group members were also new to Git, so sharing the screen during our meeting helped them to visually learn what the approval process is like on my side, and that influenced how further pull requests were made. One thing I still find challenging is how to approve parts of someone's request, not all. In that way, we could only copy-paste the text we did want to keep, reject the pull request, and I would go in and add those edits to the original. It also got tricky when stacking multiple pull requests, I do wish there was a feature where we could see side-by-side more than one pull request compared to original. Other than that, our group learned a lot together about push and pull requests straight from the command line in VS code. It was also helpful for me to learn about how to write pseudocode. 

## Other members
### Aamna
The biggest learning point for me was understanding the GitHub workflow and how forks, branches, commits, and pull requests all connect. It took some trial and error to make sure I was working in the correct branch and not overwriting anyone else’s work. I also found that reviewing multiple pull requests separately can be a little slow, especially when several people have made changes to the same file. I feel a lot more comfortable with the process now and understand much better how group coding projects are managed through GitHub.

### Selin
Like the rest of my group members, it was my first time using GitHub, so that part of this project took the most time to conceptually understand how it works and how to use it collaboratively. After my first pull request, it became a bit easier, and I can see that throughout the course, it will hopefully start to become second nature. The most difficult part was understanding how to deal with merge conflicts, when two of us had commits with conflicting code. Walking through it together on a teams call was helpful to see it from the collaborator side, however, I think once I am a group leader on one of the upcoming projects, I will even better understand the other perspective. Ultimately, this project was a great introduction to GitHub and collaborative work, as the python program itself wasn't too difficult, it allowed us room to understand the bigger picture. 

# Generative AI Appendix
Claude was used to troubleshoot accessing GitHub from the terminal and for final syntax edits of the project01.py code. 
