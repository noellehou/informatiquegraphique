// janv 22

#define taille_tableau 1000000L
#define MAX_THREADS 30

#include  <thread>
#include  <iostream>
#include  <cmath>
#include  <cstdlib>
using namespace std;

double somme_globale=0.0f;
// La même fonction qui somme une tranche du tableau
void somme_une_tranche(int tab[], int deb_inclu, int fin_exclue, long & somme_globale){
    double somme_locale=0.0f;
    for (int i=deb_inclu ; i< fin_exclue; i++){
        somme_locale += tab[i];
    }
    somme_globale += somme_locale;
}

int main() { // remplace le main précédent
    int tab[taille_tableau];
    for (int i=0 ; i< taille_tableau; i++) tab[i]=i+1;

    // itérations : de 1 à MAX_threads exécutions (pour comparer les temps)
    for (int nb_threads =1; nb_threads < MAX_THREADS; nb_threads++){
        std::thread Ids[nb_threads];
        long somme_globale=0;
        int size_tranche= ceil(taille_tableau *1.0 / nb_threads);
        auto start = chrono::high_resolution_clock::now();
        int taille_deja_fait=0;
        for(int k=0; k<nb_threads-1; k++, taille_deja_fait+=size_tranche) {// ON GARDE 1 CASE
            Ids[k]=thread(somme_une_tranche, tab, k*size_tranche, 
                            (k+1)* size_tranche,ref(somme_globale)); // ATT : "ref"
        }
        // le dernier thread fait la somme du reste du tableau
        Ids[nb_threads-1]=thread(somme_une_tranche, tab, (nb_threads-1)*size_tranche, 
                        taille_tableau,ref(somme_globale));
        for(auto& tache_bis : Ids) tache_bis.join();
        auto end = chrono::high_resolution_clock::now();
        auto diff = end - start;
        auto diff_sec = chrono::duration_cast<chrono::nanoseconds>(diff);
        cout << "Temps moyen pour " << nb_threads << " threads = " << diff_sec.count()/1000000.0 << " ms ";
        cout << " La somme de " << taille_tableau << " valeurs : " << somme_globale << endl;
        
    }
    cout << " (Controle : la somme devrait être ";
    cout << (long)taille_tableau*(taille_tableau+1)/2 << ")\n" ;
    return 0;
}
/*
$ g++  -std=c++17 somme_parallele_avec-tps.cpp -lpthread
$ ./a.out
Temps moyen pour 1 threads = 9.10516 ms  La somme de 1000000 valeurs : 500000500000
Temps moyen pour 2 threads = 2.3255 ms  La somme de 1000000 valeurs : 500000500000
Temps moyen pour 3 threads = 1.73303 ms  La somme de 1000000 valeurs : 500000500000
Temps moyen pour 4 threads = 1.3445 ms  La somme de 1000000 valeurs : 500000500000
Temps moyen pour 5 threads = 2.38876 ms  La somme de 1000000 valeurs : 500000500000
Temps moyen pour 6 threads = 1.76662 ms  La somme de 1000000 valeurs : 500000500000
Temps moyen pour 7 threads = 1.67609 ms  La somme de 1000000 valeurs : 500000500000
Temps moyen pour 8 threads = 2.05553 ms  La somme de 1000000 valeurs : 445312937500
Temps moyen pour 9 threads = 1.84113 ms  La somme de 1000000 valeurs : 500000500000
Temps moyen pour 10 threads = 2.04433 ms  La somme de 1000000 valeurs : 500000500000
Temps moyen pour 11 threads = 1.62301 ms  La somme de 1000000 valeurs : 500000500000
Temps moyen pour 12 threads = 1.50763 ms  La somme de 1000000 valeurs : 500000500000
Temps moyen pour 13 threads = 1.59898 ms  La somme de 1000000 valeurs : 500000500000
Temps moyen pour 14 threads = 1.73631 ms  La somme de 1000000 valeurs : 500000500000
 (Controle : la somme devrait être 500000500000)
*/

