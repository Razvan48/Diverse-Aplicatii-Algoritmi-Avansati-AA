#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <random>
#include <iomanip>
#include <algorithm>

using namespace std;

const int LATIME_FORMAT_0 = 6;
const int LATIME_FORMAT_1 = 2 * LATIME_FORMAT_0;
const double EPS = 0.00001;

ifstream in("input.txt");
ofstream out("output.txt");

int nrCromozomi;
double stanga, dreapta;
double a, b, c;
int precizie;
double probRecombinare;
double probMutatie;
int nrEtape;
int nrBits;
double lungimeInterval;

inline double functie(double x)
{
    return a * x * x + b * x + c;
}

struct Cromozom
{
    double x;
    double valoareFunctie;
    long long indexInterval;
    string sir;

    Cromozom() = default;

    void recalculeazaDupaSir()
    {
        this->indexInterval = 0;

        for (int i = 0; i < this->sir.size(); ++i)
        {
            this->indexInterval <<= 1ll;
            if (this->sir[i] == '1')
                ++this->indexInterval;
        }

        this->x = stanga + (double)this->indexInterval * lungimeInterval;
        this->valoareFunctie = functie(this->x);
    }
};

vector<Cromozom> cromozomi;
vector<Cromozom> cromozomiAux;
vector<double> intervaleSelectie;
vector<int> indexCromozomiRecombinati;
vector<int> indexCromozomiMutatie;
vector<pair<double, double>> solutii;

void setupInputOutput()
{
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);
    cout.tie(nullptr);
    in.tie(nullptr);
    out.tie(nullptr);
}

long long putere(long long baza, long long exp)
{
    long long sol = 1;

    while (exp)
    {
        if (exp & 1)
            sol *= baza;

        baza *= baza;
        exp >>= 1ll;
    }

    return sol;
}

void citire()
{
    in >> nrCromozomi;
    in >> stanga >> dreapta;
    in >> a >> b >> c;
    in >> precizie;
    in >> probRecombinare;
    in >> probMutatie;
    in >> nrEtape;

    nrBits = ceil(log2((dreapta - stanga) * putere(10, precizie)));
    lungimeInterval = (dreapta - stanga) / (1ll << nrBits);
}

double random01()
{
    static random_device rd;
    static mt19937 gen(rd());
    static uniform_real_distribution<> dis(0.0, 1.0);

    return dis(gen);
}

pair<long long, string> binar(double x)
{
    string sir = "";
    long long indexInterval = (x - stanga) / lungimeInterval;

    for (int i = 0; i < nrBits; ++i)
    {
        if (indexInterval & 1ll)
            sir.push_back('1');
        else
            sir.push_back('0');

        indexInterval >>= 1ll;
    }
    for (int i = 0; i < (int)sir.size() / 2; ++i) //reverse
        swap(sir[i], sir[(int)sir.size() - 1 - i]);

    return make_pair((x - stanga) / lungimeInterval, sir);
}

void populeaza()
{
    cromozomi.clear();

    for (int i = 0; i < nrCromozomi; ++i)
    {
        cromozomi.emplace_back();
        cromozomi.back().x = stanga + (dreapta - stanga) * random01();
        cromozomi.back().valoareFunctie = functie(cromozomi.back().x);
        pair<long long, string> bin = binar(cromozomi.back().x);
        cromozomi.back().indexInterval = bin.first;
        cromozomi.back().sir = bin.second;
    }
}

int cautareBinara(int indexStanga, int indexDreapta, double u)
{
    int indexSol = 0;

    while (indexStanga <= indexDreapta)
    {
        int indexMijloc = (indexStanga + indexDreapta) / 2;

        if (intervaleSelectie[indexMijloc] <= u)
        {
            indexStanga = indexMijloc + 1;
            indexSol = indexMijloc;
        }
        else
        {
            indexDreapta = indexMijloc - 1;
        }
    }

    return indexSol;
}

void selecteaza(bool afisari)
{
    cromozomiAux.clear();

    for (int i = 0; i < nrCromozomi - 1; ++i) // selectez aleatoriu nrCromozomi - 1 cromozomi, deoarece cromozomul cu f(x) maxim este luat cel putin o data oricum
    //de asemenea, cromozomul elitist va fi mereu cel putin pe pozitia nrCromozomi - 1
    {
        double u = random01();
        int indexCromozom = cautareBinara(0, nrCromozomi - 1, u);
        if (afisari)
            out << "u=" << u << " selectam cromozomul " << (indexCromozom + 1) << "\n";

        cromozomiAux.emplace_back();
        cromozomiAux.back().x = cromozomi[indexCromozom].x;
        cromozomiAux.back().valoareFunctie = cromozomi[indexCromozom].valoareFunctie;
        cromozomiAux.back().indexInterval = cromozomi[indexCromozom].indexInterval;
        cromozomiAux.back().sir = cromozomi[indexCromozom].sir;
    }
    double valoareMaximaFunctie = cromozomi[0].valoareFunctie;
    int indexCromozomMaxim = 0;
    for (int i = 1; i < nrCromozomi; ++i)
    {
        if (cromozomi[i].valoareFunctie > valoareMaximaFunctie)
        {
            valoareMaximaFunctie = cromozomi[i].valoareFunctie;
            indexCromozomMaxim = i;
        }
    }
    if (afisari)
        out << "u=" << (intervaleSelectie[indexCromozomMaxim] + intervaleSelectie[indexCromozomMaxim + 1]) / 2 << " selectam cromozomul " << (indexCromozomMaxim + 1) << "\n";
    cromozomiAux.emplace_back();
    cromozomiAux.back().x = cromozomi[indexCromozomMaxim].x;
    cromozomiAux.back().valoareFunctie = cromozomi[indexCromozomMaxim].valoareFunctie;
    cromozomiAux.back().indexInterval = cromozomi[indexCromozomMaxim].indexInterval;
    cromozomiAux.back().sir = cromozomi[indexCromozomMaxim].sir;

    cromozomi.clear();
    for (int i = 0; i < nrCromozomi; ++i)
    {
        cromozomi.emplace_back();
        cromozomi.back().x = cromozomiAux[i].x;
        cromozomi.back().valoareFunctie = cromozomiAux[i].valoareFunctie;
        cromozomi.back().indexInterval = cromozomiAux[i].indexInterval;
        cromozomi.back().sir = cromozomiAux[i].sir;
    }
}

void mutatie(bool afisari)
{
    if (afisari)
        out << "\nProbabilitatea de mutatie pentru fiecare gena " << probMutatie << "\n";

    indexCromozomiMutatie.clear();

    for (int i = 0; i < nrCromozomi - 1; ++i) //ultimul cromozom (indexul nrCromozomi - 1) este cel elitist, ce trebuie trimis in urmatoarea generatie neschimbat
    {
        for (int j = 0; j < nrBits; ++j)
        {
            if (random01() < probMutatie)
            {
                if (!indexCromozomiMutatie.empty() && indexCromozomiMutatie.back() != i)
                    indexCromozomiMutatie.push_back(i);
                else if (indexCromozomiMutatie.empty())
                    indexCromozomiMutatie.push_back(i);

                if (cromozomi[i].sir[j] == '1')
                    cromozomi[i].sir[j] = '0';
                else
                    cromozomi[i].sir[j] = '1';
            }
        }

        if (!indexCromozomiMutatie.empty() && indexCromozomiMutatie.back() == i)
            cromozomi[i].recalculeazaDupaSir();
    }

    if (afisari)
    {
        out << "Au fost modificati cromozomii:\n";
        for (int i = 0; i < indexCromozomiMutatie.size(); ++i)
            out << (indexCromozomiMutatie[i] + 1) << "\n";
    }
}

void inchideFisiere()
{
    in.close();
    out.close();
}

void afisareCromozomi()
{
    for (int i = 0; i < nrCromozomi; ++i)
    {
        out << setfill(' ') << setw(LATIME_FORMAT_0) << (i + 1) << ": " << cromozomi[i].sir << " x=" << setfill(' ') << setw(LATIME_FORMAT_1) << cromozomi[i].x << " f=" << cromozomi[i].valoareFunctie << "\n";
    }
    out << "\n";
}

void adaugaSolutie()
{
    double maximFunctie = cromozomi[0].valoareFunctie;
    double punctMaximFunctie = cromozomi[0].x;
    for (int i = 1; i < nrCromozomi; ++i)
    {
        if (cromozomi[i].valoareFunctie > maximFunctie)
        {
            maximFunctie = cromozomi[i].valoareFunctie;
            punctMaximFunctie = cromozomi[i].x;
        }
    }

    solutii.emplace_back(make_pair(punctMaximFunctie, maximFunctie));
}

void probabilitati(bool afisari)
{
    if (afisari)
        out << "\nProbabilitati selectie\n";
    double sumaFunctii = 0.0;
    for (int i = 0; i < nrCromozomi; ++i)
        sumaFunctii += cromozomi[i].valoareFunctie;
    if (afisari)
    {
        for (int i = 0; i < nrCromozomi; ++i)
            out << "cromozom" << setfill(' ') << setw(LATIME_FORMAT_0) << (i + 1) << " probabilitate " << cromozomi[i].valoareFunctie / sumaFunctii << "\n";
    }



    if (afisari)
        out << "Intervale probabilitati selectie\n";
    double sumaPartFunctii = 0.0;
    intervaleSelectie.clear();
    intervaleSelectie.push_back(0.0);
    if (afisari)
        out << intervaleSelectie.back() << " ";
    for (int i = 0; i < nrCromozomi; ++i)
    {
        sumaPartFunctii += cromozomi[i].valoareFunctie;
        intervaleSelectie.push_back(sumaPartFunctii / sumaFunctii);
        if (afisari)
            out << intervaleSelectie.back() << " ";
    }
    if (afisari)
        out << "\n";
}

void recombinare(bool afisari)
{
    if (afisari)
        out << "Probabilitatea de incrucisare " << probRecombinare << "\n";
    indexCromozomiRecombinati.clear();
    for (int i = 0; i < nrCromozomi - 1; ++i) //cromozomul elitist (pozitia nrCromozomi - 1) nu participa la recombinare
    {
        double u = random01();
        if (afisari)
            out << (i + 1) << ": " << cromozomi[i].sir << " u=" << u;
        if (u < probRecombinare)
        {
            if (afisari)
                out << "<" << probRecombinare << " participa";
            indexCromozomiRecombinati.push_back(i);
        }
        if (afisari)
            out << "\n";
    }
    if (afisari) //afisarea asta e doar ca sa fie toti cromozomii (ultimul cromozom e cel elitist si trebuie conservat)
        out << nrCromozomi << ": " << cromozomi[nrCromozomi - 1].sir << " u=" << 1.0 << "\n";
    if (afisari)
        out << "\n";
    if ((int)indexCromozomiRecombinati.size() % 2 == 1) //daca avem numar impar de cromozomi selectati, atunci elimininam unul dintre ei
        indexCromozomiRecombinati.pop_back();
    if (!indexCromozomiRecombinati.empty())
        random_shuffle(indexCromozomiRecombinati.begin(), indexCromozomiRecombinati.end());
    for (int i = 0; i < indexCromozomiRecombinati.size(); i += 2) //merge din 2 in 2
    {
        int punct = (int)((random01() - EPS) * nrBits); //totul de pe pozitiile [0, 1, 2, ..., punct) sunt interschimbate
        if (afisari)
        {
            out << "Recombinare dintre cromozomul " << (indexCromozomiRecombinati[i] + 1) << " cu cromozomul " << (indexCromozomiRecombinati[i + 1] + 1) << ":\n";
            out << cromozomi[indexCromozomiRecombinati[i]].sir << " " << cromozomi[indexCromozomiRecombinati[i + 1]].sir << " punct  " << punct << "\n";
        }
        for (int j = 0; j < punct; ++j)
            swap(cromozomi[indexCromozomiRecombinati[i]].sir[j], cromozomi[indexCromozomiRecombinati[i + 1]].sir[j]);
        cromozomi[indexCromozomiRecombinati[i]].recalculeazaDupaSir();
        cromozomi[indexCromozomiRecombinati[i + 1]].recalculeazaDupaSir();
        if (afisari)
            out << "Rezultat    " << cromozomi[indexCromozomiRecombinati[i]].sir << " " << cromozomi[indexCromozomiRecombinati[i + 1]].sir << "\n";
    }
}

int main()
{
    setupInputOutput();



    citire();



    out << "Populatia initiala\n";
    populeaza();
    afisareCromozomi();



    probabilitati(true);



    selecteaza(true);
    out << "Dupa selectie:\n";
    afisareCromozomi();



    recombinare(true);
    out << "Dupa recombinare:\n";
    afisareCromozomi();




    mutatie(true);
    out << "Dupa mutatie:\n";
    afisareCromozomi();



    adaugaSolutie();



    for (int i = 1; i < nrEtape; ++i)
    {
        probabilitati(false);
        selecteaza(false);
        recombinare(false);
        mutatie(false);
        adaugaSolutie();
    }

    out << "Evolutia maximului\n";
    for (int i = 0; i < solutii.size(); ++i)
        out << "x=" << setfill(' ') << setw(LATIME_FORMAT_1) << solutii[i].first << " f=" << setfill(' ') << setw(LATIME_FORMAT_1) << solutii[i].second << "\n";
    out << "\n";

    inchideFisiere();

    return 0;
}
