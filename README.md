# Powder_Metal_CREEP
Duva&amp;Crow (1994) pressure dependent creep model implemented with Abaqus CREEP subroutine<br><br>
Written on Intel® Fortran Version 19.0.5 (Abaqus 2023)

# AUTHOR
Name: Youngbin LIM
<br>E-mail: lyb0684@naver.com<br>

# Run command
abaqus job=model_name user=Powder_Creep.for int cpus=(number of cpus)
<br>

# Single Element Test
Single element test was performed for in-plane compression with plain strain condition to validate the CREEP sbroutine. The relative density history was compared with reference case (Introduction to Computational Plasticity, Dunne and Petrinic, 2005). Material parameters for the validation model is applied in the subroutine and input files<br><br>
![Single_Element](https://github.com/YB-LIM/Powder_Metal_CREEP/assets/105615106/312ab1fb-3758-4def-900a-0b6c2f22cd84)

# Indentation Simulation
Countours from Indentation simulation 
![Abstract](https://github.com/YB-LIM/Powder_Metal_CREEP/assets/105615106/dd6e4ed1-e905-456a-ad1f-304b81bd759d)

# Powder Compaction Forming of a Gear
Countours from compaction forming of a gear simulation
![Honeycam 2023-06-28 15-22-43](https://github.com/YB-LIM/Powder_Metal_CREEP/assets/105615106/5f597651-976a-4081-80b4-3fa9a055a42e)
