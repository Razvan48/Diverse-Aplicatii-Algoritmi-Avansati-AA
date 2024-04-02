#include <iostream>
#include <random>

using namespace std;

double random01()
{
    static random_device rd;
    static mt19937 gen(rd());
    static uniform_real_distribution<> dis(0.0, 1.0);

    return dis(gen);
}

int main()
{
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);
    cout.tie(nullptr);

    int nrTeste;
    cin >> nrTeste;

    int nrPuncteInCerc = 0;

    for (int i = 0; i < nrTeste; ++i)
    {
        double x = random01() * 2.0 - 1.0;
        double y = random01() * 2.0 - 1.0;

        if (x * x + y * y <= 1.0)
            ++nrPuncteInCerc;
    }

    double prob = 1.0 * nrPuncteInCerc / nrTeste;

    cout << prob << '\n';
    cout << "PI = " << 4.0 * prob << '\n';

    return 0;
}
