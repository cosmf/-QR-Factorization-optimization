# QR Factorization

**Titlul proiectului:** Factorizarea QR 

**Descrierea proiectului:**

    - Proiectul pornește de la implementarea serială a factorizării QR a unei matrice pătratice folosind reflectori Householder. Factorizarea implică descompunerea unei matrice A în matricele Q (ortogonală) și R (triunghiulară superioară).

    - Vom paraleliza mai multe etape pentru a îmbunătăți performanța, folosind biblioteca OpenMP:

        1. Calculul vectorului reflector Householder – formarea vectorului și calculul normei vor fi împărțite pe fire de execuție independente. 

        2. Actualizarea matricilor R și Q – fiecare coloană/rând va fi actualizată în paralel pentru a reduce timpul de execuție.

    - Paralelizarea va îmbunătăți semnificativ eficiența, în special pentru matrici mari, reducând timpul de execuție și utilizând mai eficient resursele CPU.

**Asistentul proiectului:** Tudor-Alexandru Calafeteanu

**Membrii echipei**: Mușat Irene Mihaela, Florescu Cosmin Mihai

**Complexități:**

 - Inițializarea matricelor R și Q: O(N²)

    Parcurgem fiecare element din matricile R și Q pentru inițializare, ceea ce implică două bucle imbricate care acoperă toate cele  elemente.

 - Formarea vectorilor v și u: O(N) pentru fiecare coloană

    Pentru fiecare coloană, vectorii v și u sunt calculați prin parcurgerea elementelor rămase din acea coloană, ceea ce implică o buclă de dimensiune N.

 - Actualizarea matricelor R și Q: O(N²) pentru fiecare pas

    Actualizăm matricile R și Q pentru fiecare reflector, iar acest proces implică două bucle imbricate care parcurg rândurile și coloanele rămase, adică aproximativ  operații pentru fiecare pas.

 - Număr total de pași (reflectori): N

    Aplicăm un reflector pentru fiecare coloană, deci numărul de pași este N.

 - Complexitatea totală a algoritmului: O(N³)

    Deoarece avem N pași și fiecare pas implică  operații pentru actualizarea matricelor, complexitatea totală este .

## CPU Arhitecture
   - Arhitecture: intel
   - Cpu model: Intel(R) Core(TM) i5-8350U
   - Cores: 8
   - Cpu max frequency:  1.7GHz
   - 16GB RAM

## Graphics

   #### *Timpii pentru codul serial:*
![Description of the image](plots/plots%20images/grafic_serial.png)

   #### *Timpii pentru Acceleratia OpenMP:*
![Description of the image](plots/plots%20images/acceleration_vs_cores_openmp.png)

   #### *Timpii pentru Acceleratia MPI:*
![Description of the image](plots/plots%20images/acceleration_vs_cores_mpi.png)

   #### *Timpul de executie pe cele 3 tipuri de date:*
![Description of the image](plots/plots%20images/execution_time_vs_matrix_dimension_hybrid.png)

   #### *Timpul de executie pentru o matrice de 800 x 800 pentru 1,..,8 thread-uri si procesoare:*
![Description of the image](plots/plots%20images/QR_plot_times_cores_hybrid_small.png)

   #### *Timpul de executie pentru o matrice de 900 x 900 pentru 1,..,8 thread-uri si procesoare:*
![Description of the image](plots/plots%20images/QR_plot_times_cores_hybrid_medium.png)

   #### *Timpul de executie pentru o matrice de 1100 x 1100 pentru 1,..,8 thread-uri si procesoare:*
![Description of the image](plots/plots%20images/QR_plot_times_cores_hybrid_large.png)


## Weeklog

### Week 0
- implemented serial version for blobs dataset
- profiled the program using intel Vtune on a Windows SO

### Week 1
- updated serial version for the data
- profiling made

### Week 2
- implemented openmp version
- openmp profiling

### Week 3
- implemented mpi version
- mpi profiling
- added some graphs

### Week 4
- implemented hybrid version (mpi + openMP)
- hybrid profiling
- added some graphs

### Week 5
- final touches
- added some style
