# Aim and scope for the **primerdiffer** package

## to build several methods to design primers for different specificity test

1. Design genome-wide specific primers for both species/sub-species/divergent sequences
    - Considering the dis-similarity between species (less than 95% in most cases)
    - Just greedy design primers and make a specifiyu check should be enough
    
    
2. Design genome-wide specific for strains/closely related sequences
    - Using one as reference and compare the other strains to the reference
    - Using the VCF file to get the insertion and deletion
    - Design primers with 3' end falls into the insertion 
            - -> primer specific for reference strain
            - -> primers specific for other strain
            
3. The pair-wised design can be extended to 3 or more sequence
    - later works