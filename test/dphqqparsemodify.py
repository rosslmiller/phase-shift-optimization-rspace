import re


# Take filename and regex for a partial wave coefficient and return coefficients
def partial_coefficient_parse(filename,regexstring):
    with open(filename, "r") as file:
        all_content = file.readlines()
    
    wave_regex = re.compile(regexstring)

    # Written to get coefficient floats and the corresponding index in the list of lines
    coefficients = []
    for i in range(len(all_content)):
        if wave_regex.findall(all_content[i]) != []:
            x = wave_regex.findall( all_content[i] )[0]
            coefficients.append( ( i ,float( [j for j in filter(None,x.split(' '))][0]) ) )

    return coefficients

# Determines how many spaces to the beginning and end of coefficients coeff
def line_create( coeff , annoyingzero = 0):
    first_string = "c01  1s0"
    last_string = ".        138.039   0."
      
    if annoyingzero:
        last_string = str(1) + last_string
    else:
        last_string = str(0) + last_string

    if(isinstance(coeff,str)):
        coeff = float(coeff)
        
    pos = coeff >= 0

    print(type(coeff))      

    if(pos):  
        coeff = f"{coeff:>12.6f}"
    else:
        coeff = f"{coeff:>12.6f}"
      
    final_str = first_string + coeff + last_string 

    return final_str


## Test linecreate()
#print(linecreate(-0.3048))

## Test partial_coefficient_parse
print( partial_coefficient_parse("dphqq.d", "c01\s*1s0\s*(-*\d+\.\d*\s*)[0,1].") )
