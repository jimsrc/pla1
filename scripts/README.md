# Scripts

---
### Generate .pdf with a list of figures of different observables of the simulations:
This generates a `.key` and a `.pdf` whith a basename given as an alphanumeric hash.
This hash stores the information on the filenames of the .h5 files that were 
analyzed (associated to the IDs given below as arguments).
```bash
./juguete.py -- --IDs 202,203,204 --legend Nm_slab,Nm_2d --prefix h_
```


---
### Quick report of the .h5 files associated to a given hash
 Generate a table (in terminal) of the parameters that are common/different between one .h5 file and another:
```bash
# e.g., this gives us two tables on terminal: one for the parameters 
# that are common, and another table with the parameters that are 
# different between one .h5 file and another.
./decode_pla1.py -- --hash 704a7a766658484f --list_contents
# To supress those tables, obviate the `--list-contents` option.
```


<!--- EOF -->
