#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <map>
#include <sstream>
#include <list>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <fstream>






using namespace std;

//structures:

typedef tuple<int, int, int> Triplet;
typedef vector<vector<Triplet>> Matrice;
typedef pair<int, int> Point;
typedef pair<float, float> Pointf;
typedef vector<pair<int, int>> Listepts;
typedef vector<pair<float, float>> Listeptsf;
typedef vector<vector<int>> MatGris;

struct droite
{
  double a, b, c; // Coefficients de l'équation ax + by + c = 0
};


Matrice initMat(int nb_lignes, int nb_colonnes){ // Crée une matrice de avec nb_lignes lignes et nb_colonnes colonnes, remplie de (0,0,0).
    Matrice M; 
    for(int i = 0;i<nb_lignes;i++){
        Triplet act{0,0,0};
        std::vector<Triplet> v(nb_colonnes,act);
        M.push_back(v);
    }
    return M;
}

MatGris initMatGris(int nb_lignes, int nb_colonnes){ // Crée une matrice en niveaux de gris avec nb_lignes lignes et nb_colonnes colonnes, remplie de (0,0,0).
    MatGris M; 
    for(int i = 0;i<nb_lignes;i++){
        std::vector<int> v(nb_colonnes,0);
        M.push_back(v);
    }
    return M;
}

void afficherMat(Matrice M){ // Affiche sur c++ une matrice
    for(unsigned int i = 0;i<M.size();i++){
        for(unsigned int j = 0;j<M[i].size();j++){
            Triplet act = M[i][j];
            int a = get<0>(act);
            int b = get<1>(act);
            int c = get<2>(act);
            std::cout<<"("<<a<<","<<b<<","<<c<<") ";
        }
        std::cout<<"\n";
    }
}

void afficherMatGris(MatGris M){
    for(unsigned int i = 0;i<M.size();i++){
        for(unsigned int j = 0;j<M[i].size();j++){
            std::cout <<M[i][j]<<" ";
        }
        std::cout<<"\n";
    }
}

//fonctions pratiques pour faire des print:

void printMap(const std::map<std::pair<float, float>, int>& myMap) {
    for (const auto& pair : myMap) {
        std::cout << "(" << pair.first.first << ", " << pair.first.second << ") : " << pair.second << std::endl;
    }
}

void printVectorOfPairs(const std::vector<std::pair<int, int>>& vec) {
    for (const auto& pair : vec) {
        std::cout << "(" << pair.first << ", " << pair.second << ")" << std::endl;
    }
}

void affiche_Liste(Listepts L){ // Affiche une liste de points
    for(unsigned int i = 0;i < L.size();i++){
        std::pair<int,int> coord = L[i];
        int x = coord.first;
        int y = coord.second;
        std::cout<<"("<<x<<","<<y<<") ";
    }
}



MatGris TransfomGris(Matrice M){ // Transforme une Matrice en MatGris
    MatGris RES = initMatGris(M.size(),M[0].size());
    for(unsigned int i = 0; i < M.size();i++){
        for(unsigned int j = 0;j<M[0].size();j++){
            Triplet act = M[i][j];
            int a = get<0>(act);
            int b = get<1>(act);
            int c = get<2>(act);
            int val = floor((a+b+c)/3);
            RES[i][j] = val;
        }
    }
    return RES;
}

Matrice lireppm(std::string lienppm){ // Transforme un fichier ppm en une matrice c++
    std::ifstream f(lienppm); // On ouvre le fichier en lecture
    Matrice M;
    if (!f.is_open()){
        std::cout << "Impossible d'ouvrir le fichier en lecture !" << std::endl;
    }
    else{
        std::string ligne1, ligne2,ligne3,ligne4,ligne5,ligne6 ; // Les 6 premières lignes du fichiers ne sont pas des valeurs numériques
        int dim1,dim2;
        getline(f,ligne1); // Les deux premières lignes ne sont pas nécessaires à la créationd de la matrice
        std::getline(f,ligne2);
        f >> dim1; // On récupère le nombre de lignes de l'image
        f >> dim2; // On récupère le nombre de colonnes de l'image
        const int d1 = dim1;
        const int d2 = dim2;
        std::getline(f,ligne3); // Les quatre lignes suivantes ne sont donc plus nécessaires
        std::getline(f,ligne4);
        std::getline(f,ligne5);
        std::getline(f,ligne6);
        M = initMat(d1,d2); // On initialise une matrice avec d1 lignes et d2 colonnes
        for(int i = 0; i < d1 * d2;i++){
            std::string valR,valG,valB; // On récupère nos valeurs de Rouge, Vert et Bleu qui sont d'abord des chaînes de caractère
            std::getline(f,valR);
            std::getline(f,valG);
            std::getline(f,valB);
            int valeurR = std::stoi(valR); // On transforme ces chaînes de caractère en entier
            int valeurG = std::stoi(valG);
            int valeurB = std::stoi(valB);
            Triplet act{valeurR,valeurG,valeurB}; // Il s'agit du prochain pixel à ajouter à notre matrice
            int ligne = i / d1; // Le numéro de ligne sur lequel on est actuellement
            int colonne = i%d1; // Le numéro de colonne sur lequel on est actuellement
            M[ligne][colonne] = act;
        }
        f.close();
    }
    return M;
}

void ecriveur(Matrice a,std::string nom){ // Ecrit la matrice en entrée dans un fichier avec le nom en entrée
  int d1 = a.size();
  int d2 = a[0].size();
  string const nomFichier(nom);
    std::ofstream monFlux(nomFichier.c_str());
    if(monFlux){
      
        monFlux << "P3" << std::endl;
        monFlux << "#resolution" << std::endl;
        monFlux << d1 << " " << d2 << std::endl;
        monFlux << "#avec 255 comme valeur max" << std::endl;
        monFlux << "255" << std::endl;
        monFlux << "#debut de l'image" << std::endl;
        for (int i = 0; i < d1; ++i){
          for (int j = 0; j < d2; ++j){
            Triplet act = a[i][j];
            monFlux << get<0>(act) << std::endl;
            monFlux << get<1>(act) << std::endl;
            monFlux << get<2>(act) << std::endl;
        }}
    }
    else
    {
        std::cout << "ERREUR: Impossible d'ouvrir le fichier." << std::endl;
    }
}

void ecriveur2(MatGris a,std::string nom){ // Ecrit la matrice en niveaux de gris en entrée dans un fichier avec le nom en entrée
  int d1 = a.size();
  int d2 = a[0].size();
  string const nomFichier(nom);
    std::ofstream monFlux(nomFichier.c_str());
    if(monFlux){
      
        monFlux << "P2" << std::endl;
        monFlux << "#Voici mon image" << std::endl;
        monFlux << d1 << " " << d2 << std::endl;
        monFlux << "255" << std::endl;
        for (int i = 0; i < d1; ++i){
            for (int j = 0; j < d2; ++j){
                int act = a[i][j];
                monFlux << act<<" ";
            }
            monFlux<<std::endl;
        }
    }
    else
    {
        std::cout << "ERREUR: Impossible d'ouvrir le fichier." << std::endl;
    }
}

Listepts NonBlancs(Matrice M){ // Renvoie la liste des coordonnées de tous les pixels qui ne sont pas blancs donc (255,255,255)
    Listepts RES;
    for(unsigned int i = 0; i <M.size();i++){
        for(unsigned int j = 0; j < M[i].size();j++){ // On parcourt les différents coefficients de la matrice
            Triplet act = M[i][j];
            int a = get<0>(act);
            int b = get<1>(act);
            int c = get<2>(act);
            if(a != 255 || b != 255 || c != 255){ // Et on vérifie s'ils sont blancs ou pas
                std::pair<int,int> coord; 
                coord.first = i;
                coord.second = j;
                RES.push_back(coord); // Si ce n'est pas le cas, on l'ajoute à notre liste resultat
            }
        }
    }
    return RES;
}

Listepts NonBlancs2(MatGris M){// Pareil que précédemment mais avec des matrices en niveaux de gris
    Listepts RES;
    for(unsigned int i = 0; i < M.size();i++){
        for(unsigned int j = 0; j < M[0].size();j++){
            if(M[i][j] != 255){
                std::cout <<i<<","<<j<<std::endl;
                std::pair<unsigned int,unsigned int> coord; 
                coord.first = i;
                coord.second = j;
                RES.push_back(coord);
            }
        }
    }
    return RES;
}

Listepts Blancs(MatGris M){ // A l'inverse de NonBlancs, renvoie la liste des points blancs de l'image
    Listepts RES;
    for(unsigned int i = 0; i < M.size();i++){
        for(unsigned int j = 0; j < M[0].size();j++){
            if(M[i][j] == 255){
                std::cout <<i<<","<<j<<std::endl;
                std::pair<unsigned int,unsigned int> coord; 
                coord.first = i;
                coord.second = j;
                RES.push_back(coord);
            }
        }
    }
    return RES;
} 

Matrice changer_en_noir(Matrice M, Listepts L){ // Modifie en noir tous les points d'une image dont les coordonnées ont été renseignées dans une liste
    Matrice RES = M;
    for(unsigned int i = 0; i < L.size();i++){
        Triplet change{0,0,0};
        RES[L[i].first][L[i].second] = change;
    }
    return RES;
}

MatGris NB(MatGris M,Listepts L){ // Transforme une matrice de niveaux de gris en une matrice en noir et blanc uniquement
    MatGris RES = M;
    for(unsigned int i = 0; i < L.size();i++){
        pair<int, int> coord = L[i];
        int x = coord.first;
        int y = coord.second;
        RES[x][y] = 0;
    }
    return RES;
}


MatGris filtreSobel(MatGris M){ // Applique un filtre de Sobel sur l'image donnée en entrée
    int pixel_x;
    int pixel_y;

    float sobel_x[3][3] =
    { { -1, 0, 1 },
      { -2, 0, 2 },
      { -1, 0, 1 } };

    float sobel_y[3][3] =
    { { -1, -2, -1 },
      { 0,  0,  0 },
      { 1,  2,  1 } };
    int nbl = M.size();
    int nbc = M[0].size();
    MatGris RES = initMatGris(nbl, nbc);
    for (int x=0; x < nbc-2; x++){
        for (int y=0; y < nbl-2; y++){
            pixel_x = (sobel_x[0][0] * M[(nbc * y+ x)/nbl][(nbc * y + x)%nbl])
                    + (sobel_x[0][1] * M[(nbc* y- +  x+1)/nbl][(nbc* y +  x+1)%nbl])
                    + (sobel_x[0][2] * M[(nbc * y + (x+2))/nbl][(nbc * y + (x+2))%nbl])
                    + (sobel_x[1][0] * M[(nbc*  (y+1)    + x)/nbl][(nbc*  (y+1)    + x)%nbl])
                    + (sobel_x[1][1] * M[(nbc *  (y+1)    +  x+1)/nbl][(nbc *  (y+1)    +  x+1)%nbl])
                    + (sobel_x[1][2] * M[(nbc *  (y+1)    + (x+2))/nbl][(nbc *  (y+1)    + (x+2))%nbl])
                    + (sobel_x[2][0] * M[(nbc * (y+2) + x)/nbl][(nbc * (y+2) + x)%nbl])
                    + (sobel_x[2][1] * M[(nbc * (y+2) +  x+1)/nbl][(nbc * (y+2) +  x +1)%nbl])
                    + (sobel_x[2][2] * M[(nbc * (y+2) + (x+2))/nbl][(nbc * (y+2) + (x+2))%nbl]);

            pixel_y = (sobel_y[0][0] * M[(nbc * y+ x)/nbl][(nbc * y + x)%nbl])
                    + (sobel_y[0][1] * M[(nbc* y- +  x+1)/nbl][(nbc* y +  x+1)%nbl])
                    + (sobel_y[0][2] * M[(nbc * y + (x+2))/nbl][(nbc * y + (x+2))%nbl])
                    + (sobel_y[1][0] * M[(nbc*  (y+1)    + x)/nbl][(nbc*  (y+1)    + x)%nbl])
                    + (sobel_y[1][1] * M[(nbc *  (y+1)    +  x+1)/nbl][(nbc *  (y+1)    +  x+1)%nbl])
                    + (sobel_y[1][2] * M[(nbc *  (y+1)    + (x+2))/nbl][(nbc *  (y+1)    + (x+2))%nbl])
                    + (sobel_y[2][0] * M[(nbc * (y+2) + x)/nbl][(nbc * (y+2) + x)%nbl])
                    + (sobel_y[2][1] * M[(nbc * (y+2) +  x+1)/nbl][(nbc * (y+2) +  x +1)%nbl])
                    + (sobel_y[2][2] * M[(nbc * (y+2) + (x+2))/nbl][(nbc * (y+2) + (x+2))%nbl]);
            int val = (int)sqrt((pixel_x * pixel_x) + (pixel_y * pixel_y));
            if(val < 0){val = 0;}
            if(val > 255){val = 255;}
            std::cout<<val<<std::endl;
            RES[(nbl * y + x)/nbl + 1][(nbl * y + x)%nbl +1] = val;
        }
    }
    return RES;
}

MatGris LtoM(Listepts L,int nbl,int nbc){// Transforme une liste de points blancs en une matrice noire et blanche
  MatGris RES = initMatGris(nbl,nbc);
  for(unsigned int i = 0; i < L.size();i++){
    std::pair<int,int> act;
    act = L[i];
    RES[act.first][act.second] = 255;
  }
  return RES;
}

Matrice surligner(Matrice M, Listepts L){
  Matrice RES = M;
  for(unsigned int i = 0; i < L.size();i++){
    std::pair<int,int> act = L[i];
    Triplet bleu{0,0,255};
    RES[act.first][act.second] = bleu;
  }
  return RES;
}
//dans l'image "matrice", colorie les pixels aux coordonnées données dans "points" en noir.
void afficherdroite(Matrice& matrice, Listepts& points) {
  for (auto& point : points) {
    unsigned int x = point.first;
    unsigned int y = point.second;

    if (x < 0 || x >= matrice.size() || y < 0 || y >= matrice[0].size()) {  //on verifie que les pointts sont bien dans la matrice (utile pour les tests et le débuggage)
        cout << "Point hors de la matrice : (" << x << ", " << y << ")" << endl;
    } else {
        matrice[x][y] = make_tuple(0, 0, 0);
        ;
    }
  }
}

//calcule puis affiche le segment qui represente une droite y=mx+b sur l'image "matrice"
void afficherdroite2(Matrice& matrice,float m,float b){
  unsigned int d2 = matrice.size();
  unsigned int d1 = matrice[0].size();
  for (unsigned int j = 0; j < 10000; ++j){ //si on travaille sur une image de plus de 10 000^2 pixels, on peut augmenter; de toute façon la boucle ne sera pas parcourue plus que la taille de la matrice grace au break.
    if (m*j+b<d1 && j < d2){
      Listepts mono {{floor(m*j+b),j}}; //on crée un des points de la droite qui doit être affichées.
      afficherdroite(matrice,mono);  //la fonction qui est au dessus et qui affiche une liste de points (ici , le seul point inclu dans "mono"
    }else{
      break;
    }
  }
}





//fonctions pour la transformée de Hough:

float distance(Point pt1,Point pt2)//calcul la distance entre deux points
{
  float res;
  res=sqrt(pow((pt1.first-pt2.first),2)+pow((pt1.second-pt2.second),2));
  return res;
};

Pointf intersection(droite r, droite s)//retournes le point d'intersection de deux droites
{
  Pointf p;
  if(r.a*s.b-s.a*r.b==0){// cas des droites parallèle, on considère une intersection à l'infini
    if(-(r.a/r.b)-(r.c/r.b)<=1){
      p.first=numeric_limits<int>::max();//On fait en sorte que le point d'intersection se situe à "l'infinie" et dans la direction des droite
      p.second=-(r.a/r.b)*(p.first)-(r.c/r.b);
        }else{//On doit faire des conditions pour éviter d'avoir une coordonnée qui dépasse max<int>
      p.second=numeric_limits<int>::max();
      p.first=-(r.b/r.a)*p.second-(r.c/r.a);
    }
    return p;
  }else{//cas classique, on rentre les coordonnées du point d'intersection de deux droites non parallèles
    p.first=(s.c*r.b-r.c*s.b)/(r.a*s.b-s.a*r.b);
    p.second=-(r.a/r.b)*((s.c*r.b-r.c*s.b)/(r.a*s.b-s.a*r.b))-(r.c/r.b);
    return p;
    }
}

droite Hough(Point p)//retourne la transformée de Hough d'un point
{
  droite h;
  h.a=-p.first;
  h.b=1;
  h.c=-p.second;
  return h;
};

Point inverseHough(droite h)//retourne la transformée inverse de Hough d'une droite //On ne s'est pas servi de cette fonction mais on la garde quand même
{
  Point p;
  p.first=-h.a;
  p.second=-h.c;
  return p;
};

Pointf centregravite(Listeptsf v)//retourne le point centre de gravité d'un vecteur de points (utilisée dans sontAlignéApprox())
{
  double t,u;
  for(unsigned int i = 0; i < v.size(); i++){//les coordonnées du centre de gravité sont la moyennes des coordonnées des points du nuage de point
    t+=v[i].first;
    u+=v[i].second;
  };
  t=t/v.size();
  u=u/v.size();
  Pointf p;
  p.first=t;
  p.second=u;
  return p;
}

Listeptsf listintersec(std::vector<droite> v)//retourne la liste des point d'intersection de toutes les droites d'un vecteurs de droites
{
  Listeptsf listpt;
  std::vector<droite> newv;
  std::vector<droite> newvshort;
  newv=v;
  for(unsigned int k=0;k<v.size();++k){
    newvshort=newv;
    newvshort.pop_back();
    for(unsigned int p=0;p<v.size()-1;++p){
      Pointf pti;
      pti=intersection(newv.back(), newvshort[p]);
      listpt.push_back(pti);
    };
    newv.pop_back();
  };
  return listpt;
};

vector<droite> espaceHough(Listepts v)//retourne le vecteur des droites issue de la transformation d'un nuage de point par la transformée de Hough
{
  std::vector<droite> E;
  for(unsigned int k=0;k<v.size();++k){
    E.push_back(Hough(v[k]));
  };
  return E;
};


Listeptsf nuageHough(Listepts v)//prends un nuage de point et retourne l'intersection des droites dans l'espace de Hough asssocié
{
  std::vector<droite> espace;
  espace=espaceHough(v);
  Listeptsf inter;
  inter=listintersec(espace);
  return inter;
}

set<Pointf> vectorToSet(Listeptsf v)//transforme un vecteur de points flottants en ensemble de points flottants, on se sert de cette fonction nottament pour supprimer les doublons lorsque l'on en a besoin
{
    set<Pointf> S;
    for(unsigned int k=0;k<v.size();k++){
        S.insert(v[k]);
    }
    return S;
}

set<Point> vectorToSetI(Listepts v)//transforme un vecteur de points d'entiers en ensemble de points d'entiers
{
    set<Point> S;
    for(unsigned int k=0;k<v.size();k++){
        S.insert(v[k]);
    }
    return S;
}

Listepts setToVector(set<Point> s)//transforme un ensemble de points d'entiers en vecteur de points d'entiers 
{
    Listepts result;
    for(const auto& element : s) {
        result.push_back(element);
    }
  return result;
}

map<Pointf, int> countOccurrences(Listeptsf vec)//retourne un dictionnaire avec les points contenue dans le vecteur en clé et leur nombre d'occurence dans ce vecteur en valeur de clé
{
    map<Pointf, int> occurrenceMap;
    for(const auto& pair : vec) {
        occurrenceMap[pair]++;
    }
    return occurrenceMap;
}

bool sontAligne(Listepts v)//indique (true ou false) si les points d'un vecteur sont alignés, cette fonction à été codé pour faire des tests initiaux mais ne sert pas dans la fonction du filtre
{
  Listeptsf W;
  W=nuageHough(v);
  set<Pointf> SW;
  SW=vectorToSet(W);
  if(SW.size()==1)//on vérifie si il y a un unique point d'intersection des droites dans l'espace de Hough
  {
    return true;
  }else{
    return false;
  };
}

Listeptsf listcles(const map<Pointf, int>& dico)//retourne la liste des clés d'un dictionnaire de type {points:occurence}
{
    Listeptsf cles;
    for (const auto& pair : dico) {
        cles.push_back(pair.first);
    }
    return cles;
}



map<Pointf, int> accumulateur(Listepts v,int n)//retourne un dictionnaire avec les points contenues dans le vecteur en clé dont le nombre d'occurence est suppérieur à n
{
  Listeptsf nuage;
  nuage=nuageHough(v);
  map<Pointf, int> dico;
  dico=countOccurrences(nuage);
  auto it = dico.begin();
    while (it != dico.end()) {//on parcours le dictionnaire par indexation jusqu'à la fin
        if (it->second <= n) {
            it = dico.erase(it);//on efface les points d'intersection de l'espace de Hough qui ne comptent pas assez d'intersections
        } else {
            ++it;//ou on passe à l'indexe suivant si la condition est respectée
        }
    }
  return dico;
}



bool intersect(droite d,Pointf point)//indique (true or false) si une droite intersecte un point 
{
  if(d.a*point.first+d.b*point.second+d.c==0){
    return true;
  }else{
    return false;
  }
}

Listepts filtreHoughpts(Listepts ptimage, int seuil)//retoune le vecteur de point fesant partie de la même droite
{
  Listepts sortiefiltre;
  map<Pointf, int> dico;
  dico=accumulateur(ptimage,seuil);
  printMap(dico);//visualisation du dictionnaire
  Listeptsf survivants;
  survivants=listcles(dico);//les points aux coordonnées des différentes accumulation qui ont été gardé
  vector<droite> lstdroite;
  lstdroite=espaceHough(ptimage);//les droites correspondantes à la transformée de Hough de chaques points
  for(unsigned int k=0;k<lstdroite.size();++k){
    for(unsigned int p=0;p<survivants.size();++p){
      if(intersect(lstdroite[k],survivants[p])){//on retrouve les points de l'image d'origine dont les droites dans l'espace de Hough intersecte les points d'intersections restantes
        sortiefiltre.push_back(ptimage[k]);//et on les concatène dans l'image de sortie
      }
    }
  }
  set<Point> Set; //supression des doublons dans le vecteur
  Set=vectorToSetI(sortiefiltre);
  Listepts sortiefiltreSansDoublon;
  sortiefiltreSansDoublon=setToVector(Set);
  return sortiefiltreSansDoublon;
};




//On a ici une fonction qui vérifie si des points sont en peu prêt alignés, qui fonctionne sur le principe que si c'est le cas alors le nuage des intersections des droites dans l'espace de Hough est inclue dans un petit disque de rayon epsilon
//à priori elle ne marche pas, surement à cause des propriétés de divergence de la transformée en coordonnée carthésienne, donc le principe est surement exploitable pour la transformée en coordonnée polaire

bool sontAligneApprox(Listepts v,float eps)//vérifie si une liste de points sont approximativement alignées
{
  Listeptsf W;
  W=nuageHough(v);
  Pointf centregrav;
  centregrav=centregravite(W);
  for(unsigned int k=0;k<v.size();++k){
    std::cout<<distance(W[k],centregrav)<<std::endl;
    if(distance(W[k],centregrav)>eps){
      return false;
    };
  };
  return true;
}



//prends aléatoirement 2 points non nuls de la matrice, calcule les coefs de la droite associée, les stocke avec un compteur d'apparition, et l'affiche sur la matrice si le seuil en entrée est atteint
Matrice CouplesAleatoires(Matrice matrice,int  nombredecouples = 10000,int  seuil = 1000){  //j'ai mis 1000 pour des exemples triviaux; pour des vraies images il faut le réduire!
    std::map<std::pair<float, float>, int> nblignes; // Pour stocker le nombre d'occurrences de chaque couple (m,p)
        
    //creation des droites (autant de droites que le second argument de la fonction "nombredecouples"
    for(int i = 0; i < nombredecouples; i++) {
      Point p1;                                                //premier point
      unsigned int k = rand() % NonBlancs(matrice).size();  
      p1.first = NonBlancs(matrice)[k].second;
      p1.second = NonBlancs(matrice)[k].first;  
      Point p2;                                                //deuxieme point
      unsigned int k2 = rand() % NonBlancs(matrice).size();
      p2.first = NonBlancs(matrice)[k2].second;
      p2.second = NonBlancs(matrice)[k2].first;       
      if (p2.first != p1.first){                          //si la droite est verticale, ça ne marche pas; donc on évite ce cas sans le traiter
        float m = (float) (p2.second - p1.second) / (p2.first - p1.first);  //on calcule les coordonnées de la droite associée aux deux points
        float p = p1.second - m * p1.first;   
        nblignes[{m, p}]=nblignes[{m, p}]+1;  //on indente le dictionnaire pour la droite qu'on vient de trouver
      }
    }
  
    //affichage ssi nbligne(droite) >= seuil
    for(auto& pair : nblignes) {
        if(pair.second >= seuil) {
            afficherdroite2(matrice,pair.first.first,pair.first.second); // remarque: la variable "pair" est de la forme {{m,p},n}
            std::cout << "Droite : y = " << pair.first.first << "x + " << pair.first.second << std::endl; 
      }
    }
    return matrice;
}








//tests

int main(){
  //test intersection:
  /*
  point p;
  droite r;
  droite s;
  r.a=1;
  r.b=2;
  r.c=3;
  s.a=4;
  s.b=5;
  s.c=6;
  p=r.intersection(r,s);
  std::cout<<p.first<<p.second<<std::endl;
  
  //test centre de gravité
  
  point A;
  point B;
  A.first=1;
  A.second=2;
  B.first=3;
  B.second=4;
  
  std::vector<point> w;
  w.push_back(A);
  w.push_back(B);
  point pt;
  pt=centregravite(w);
  std::cout<<pt.first<<std::endl;
  */
  
  //On introduit un vecteur de point, dont 4 sont alignés, 3 de même, et 1 n'est aligné avec aucun autre
  Point C;
  Point D;
  Point E;
  Point F;
  Point G;
  Point H;
  Point I;
  C.first=2;
  C.second=1;
  D.first=4;
  D.second=3;
  E.first=6;
  E.second=5;
  F.first=1;
  F.second=2;
  G.first=8;
  G.second=7;
  H.first=3;
  H.second=0;
  I.first=6;
  I.second=1;
  
  Listepts s;
  s.push_back(C);
  s.push_back(D);
  s.push_back(E);
  s.push_back(F);
  s.push_back(G);
  s.push_back(H);
  s.push_back(I);
  Listepts f;
  f=filtreHoughpts(s,10); //on met un seuil élevé pour retenir seulement un segment
  printVectorOfPairs(f);
  /*
  //test d'intersection
  droite A;
  droite B;
  droite C;
  A.a=1;
  A.b=1;
  A.c=0;
  B.a=1;
  B.b=2;
  B.c=1;
  C.a=4;
  C.b=3;
  C.c=3;
  std::vector<droite> v;
  v.push_back(A);
  v.push_back(B);
  v.push_back(C);
  Pointf lp;
  lp=intersection(B,C);
  std::cout<<lp.first<<std::endl;
  */
    Matrice  M;
    M = initMat(7,5); // Crée une matrice remplie de 0 de 7 lignes et 5 colonnes
    afficherMat(M); // Affiche M
    Matrice image;
    image = lireppm("Projet/m1projetcpp2.ppm"); // Stocke notre image d'un M1 rouge
    afficherMat(image); // En affiche la matrice
    Matrice image2;
    image2 = lireppm("Projet/imgSimple.ppm");
    afficherMat(image2);

  // Test de afficherdroite
    Matrice matrice = {{{1, 2, 3}, {4, 5, 6}}, {{7, 8, 9}, {10, 11, 12}}};
    afficherMat(matrice);
    Listepts points = {{0, 0}, {1, 1}};
    afficherdroite(matrice, points);
    afficherMat(matrice);

    Listepts L;
    L = NonBlancs(image);
    affiche_Liste(L);
    Matrice image4 = changer_en_noir(image,L);  // Change en noir les points renseignés dans la liste
    afficherMat(image4);
    ecriveur(image4,"Projet/ImageNoirBlanc.ppm"); 
    MatGris image5 = TransfomGris(image);
    afficherMatGris(image5);
    ecriveur2(image5,"Projet/ImageGris.ppm"); // Trasforme notre image en niveaux de gris
    MatGris image6 = filtreSobel(image5); // Applique le filtre de sobel sur cette dernière image
    afficherMatGris(image6);
    ecriveur2(image6,"Projet/ImageSobel.ppm");
    Listepts L2 = NonBlancs2(image6);
    affiche_Liste(L2);
    MatGris image7 = NB(image6,L2); // Transforme image6 en noir et blanc 
    afficherMatGris(image7);
    ecriveur2(image7,"Projet/SobelNB.ppm");

    Listepts L3 = Blancs(image7);
    affiche_Liste(L3);
    Listepts L4 = filtreHoughpts(L3,75); // Applique le filtre de Hough sur notre dernière image
    affiche_Liste(L4);
    MatGris image8 = LtoM(L4,30,30);
    ecriveur2(image8,"Projet/HoughNB.ppm");
    Matrice image9(surligner(image,L4));
    ecriveur(image9,"Projet/surligner.ppm");

    Matrice image10 = CouplesAleatoires(image,10000,20);
    afficherMat(image10);
    ecriveur(image10,"Projet/alea.ppm");


  return 0;
}
