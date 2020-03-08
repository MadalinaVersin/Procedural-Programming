#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef struct  pixel {

   unsigned char b;
   unsigned char g;
   unsigned char r;
} pixel;

typedef struct detectie {

    unsigned int linie;
    unsigned int coloana;
    double corelatie;
    pixel culoare;
} detectie;

unsigned int * xorshift32(unsigned int seed, unsigned int n) {

	unsigned int r, i;
    unsigned int *v;
    v = (unsigned int * )malloc(n * sizeof(unsigned int));
    r = seed;
    v[0] = seed;
    for (i = 1; i < n; i++) {
        r ^= r << 13;
        r ^= r >> 17;
        r ^= r << 5;
        v[i] = r;
    }

	return v;
}
unsigned int *permutare(unsigned int *array, unsigned int n) {

    unsigned int* sigma;
    unsigned int i;
    sigma = (unsigned int *)malloc(n * sizeof(unsigned int));

    for (i = 0; i < n; i++) {
        sigma[i]=i;
    }

    unsigned int p = 0;
    unsigned int j, tmp;

    for (i = n - 1; i > 0; i--) {
        j = array[p]%(i + 1);
        tmp = sigma[i];
        sigma[i] = sigma[j];
        sigma[j] = tmp;
        p = p + 1;
   }

   return sigma;
}

pixel *liniarizare(char* nume_fisier_sursa, unsigned int *n, unsigned int*m) {
    pixel *img_lin;
    unsigned int dim_img, latime_img, inaltime_img;
    FILE *fin = fopen (nume_fisier_sursa, "rb");
    if (fin == NULL) {

        printf("Nu a putut fi deschisa imaginea sursa din care citesc!");
        return ;

    }

    fseek(fin, 2, SEEK_SET);
    fread(&dim_img, sizeof(unsigned int), 1, fin);
    fseek(fin, 18, SEEK_SET);
    fread(&latime_img, sizeof(unsigned int), 1, fin);
    fread(&inaltime_img, sizeof(unsigned int), 1, fin);
    *n = inaltime_img;
    *m = latime_img;
    int padding;

    if (latime_img % 4 != 0) {
         padding=4-(3*latime_img)%4;
    }
    else {
         padding=0;
    }

    img_lin = (pixel*)malloc(latime_img * inaltime_img * (sizeof(pixel)));
    fseek(fin, 54, SEEK_SET);
    unsigned int i, j, k, l;
    unsigned char c;
    k = 0;
    for (i = 0; i < inaltime_img; i++) {
        for (j = 0; j < latime_img; j++) {
                    for(l = 0; l < 3; l++) {
                        fread(&c, 1, 1, fin);

                        if ( l == 0){
                            img_lin[k].b=c;
                        }
                        else if ( l==1 ) {

                                img_lin[k].g=c;
                                }
                                else if (l == 2) {

                                        img_lin[k].r=c;
                                }
                    }
                    k = k + 1;
        }
        fseek (fin, padding, SEEK_CUR);
    }

   fclose (fin);
   return img_lin;

}
pixel *liniarizareinversa(pixel* img, unsigned int h, unsigned int w) {

    int i1, i2, j;
    pixel aux;
    i1 = 0;
    i2 = w * h - w;

    while (i1 < i2) {

        for ( j = 0; j < w; j++) {

                aux.r = img[i1+j].r;
                aux.g = img[i1+j].g;
                aux.b = img[i1+j].b;
                img[i1+j].r = img[i2+j].r;
                img[i1+j].g = img[i2+j].g;
                img[i1+j].b = img[i2+j].b;
                img[i2+j].r = aux.r;
                img[i2+j].g = aux.g;
                img[i2+j].b = aux.b;
        }
        i1 = i1 + w;
        i2 = i2 - w;

    }
  return img;

}

void salveaza_extern(pixel* img_lin, char* nume_fisier_sursa, char*nume_fisier_destinatie) {

    unsigned int dim_img, inaltime_img, latime_img, i;
    FILE *fout = fopen(nume_fisier_destinatie, "wb");
    FILE *fin = fopen(nume_fisier_sursa, "rb");

    if((fout==NULL)||(fin==NULL)) {

        printf("Fisierul in care dorim sa salvam extern nu poate fi deschis!\n");
        return;
    }

    fseek(fin, 18, SEEK_SET);
    fread(&latime_img, sizeof(unsigned int), 1, fin);
    fread(&inaltime_img, sizeof(unsigned int), 1, fin);
    fseek(fin, 0, SEEK_SET);

    unsigned char c;
    while (fread(&c, 1, 1, fin) == 1) {

            fwrite(&c, 1, 1, fout);
		    fflush(fout);

	}
	int padding;
    if (latime_img % 4 != 0) {

         padding = 4 - (3 * latime_img) % 4;
    }
    else {
         padding=0;
    }

	fseek(fout, 54, SEEK_SET);

    for (i=0; i < inaltime_img * latime_img; i++) {

        fwrite(&img_lin[i].b, sizeof(unsigned char), 1, fout);
        fflush(fout);
        fwrite(&img_lin[i].g, sizeof(unsigned char), 1, fout);
        fflush(fout);
        fwrite(&img_lin[i].r, sizeof(unsigned char), 1, fout);
        fflush(fout);

        if ((i + 1) % latime_img == 0) {

                fseek(fout,padding,SEEK_CUR);
        }
     }

     fclose(fout);
     fclose(fin);
}

void criptare(char *fisier_sursa, char*fisier_destinatie, char*fisier_cheie) {
    unsigned int *r;
    pixel *img, *imgnou, *img_per, *img_inv;
    unsigned int w, h, sv;
    FILE *fin = fopen(fisier_sursa, "rb");
    FILE *fout = fopen(fisier_destinatie, "wb");
    FILE *f = fopen(fisier_cheie, "r");

    if (fin == NULL || fout == NULL || f == NULL) {

         printf("Unul din fisiere nu poate fi deschis!\n");
         return;
     }

    unsigned int seed, i;

    fscanf(f, "%u %u", &seed, &sv);

    img = liniarizare(fisier_sursa, &h, &w);
    img_inv = liniarizareinversa(img, h, w);
    r = xorshift32(seed, 2 * w * h);
    unsigned int *sigma;
    sigma = permutare(r + 1, h * w);
    imgnou = (pixel*)malloc(h * w * sizeof(pixel));

    for (i = 0; i < h* w; i++) {
            imgnou[sigma[i]].r = img_inv[i].r;
            imgnou[sigma[i]].g = img_inv[i].g;
            imgnou[sigma[i]].b = img_inv[i].b;
    }

    for (i = 0; i < h * w; i++) {
            printf(" %u %u %u ",imgnou[i].r,imgnou[i].g,imgnou[i].b);
    }

    pixel *imgcript;
    unsigned char *psv;
    psv = &sv;
    imgcript = (pixel*)malloc(h * w * sizeof(pixel));
    for (i = 0; i < h * w; i++) {

            unsigned char *pr;
            pr = &r[i + h * w];
            if (i == 0) {
                imgcript[i].r = (*(psv+2)) ^ imgnou[i].r ^ (*(pr + 2));
                imgcript[i].g = (*(psv + 1)) ^ imgnou[i].g ^ (*(pr + 1));
                imgcript[i].b = (*(psv + 0)) ^ imgnou[i].b ^ (*(pr + 0));
            }
           else {

                imgcript[i].r = imgcript[i - 1].r ^ imgnou[i].r ^ (*(pr + 2));
                imgcript[i].g = imgcript[i - 1].g ^ imgnou[i].g ^ (*(pr + 1));
                imgcript[i].b = imgcript[i - 1].b ^ imgnou[i].b ^ (*(pr + 0));
           }
     }

    liniarizareinversa(imgcript, h, w);
    salveaza_extern(imgcript, fisier_sursa, fisier_destinatie);
    fclose(fin);
    fclose(fout);
    fclose(f);
}

void decriptare(char *fisierimg, char*fisierimgcriptata, char*cheiasecreta) {

    FILE *fimg = fopen(fisierimg, "wb");
    FILE *fcript = fopen(fisierimgcriptata, "rb");
    FILE *f = fopen(cheiasecreta, "r");

    if ((fimg == NULL) || (fcript == NULL) || (f == NULL)) {
        printf("Unul din fisiere nu poate fi deschis!\n");
        return;
    }

    unsigned int h, w, i;
    fseek(fcript, 18, SEEK_SET);
    fread(&w, sizeof(unsigned int), 1, fcript);
    fread(&h, sizeof(unsigned int), 1, fcript);

    unsigned int seed, sv;
    fscanf(f, "%u %u", &seed, &sv);
    unsigned int *r, *sigma, *invsigma;
    pixel *imgcript, *imgint;

    r = xorshift32(seed, 2 * w * h);
    sigma = permutare(r + 1, h * w);

    invsigma = (unsigned int*)malloc(h * w * (sizeof(unsigned int)));

    for (i = 0; i < h * w; i++) {
            invsigma[sigma[i]]=i;
    }

    imgcript = (pixel*)malloc(h * w * sizeof(pixel));
    imgcript = liniarizare(fisierimgcriptata, &h, &w);
    liniarizareinversa(imgcript, h, w);

    unsigned char *psv, *pr;
    imgint = (pixel*)malloc(h * w * sizeof(pixel));
    psv = &sv;
    for (i = 0; i < h * w; i++) {

        pr = &r[i + h * w];

        if (i == 0) {

              imgint[i].r = (*(psv + 2)) ^ imgcript[i].r ^ (*(pr + 2));
              imgint[i].g = (*(psv + 1)) ^ imgcript[i].g ^ (*(pr + 1));
              imgint[i].b = (*(psv + 0)) ^ imgcript[i].b ^ (*(pr + 0));
        }
        else {
                imgint[i].r = imgcript[i - 1].r ^ imgcript[i].r ^ (*(pr + 2));
                imgint[i].g = imgcript[i - 1].g ^ imgcript[i].g ^ (*(pr + 1));
                imgint[i].b = imgcript[i - 1].b ^ imgcript[i].b ^ (*(pr + 0));
        }
    }

    pixel *d;
    d = (pixel*)malloc(h * w * sizeof(pixel));
    for (i = 0; i < h * w; i++) {
         d[invsigma[i]].r = imgint[i].r;
         d[invsigma[i]].g = imgint[i].g;
         d[invsigma[i]].b = imgint[i].b;
     }
    liniarizareinversa(d, h, w);
    salveaza_extern(d, fisierimgcriptata, fisierimg);
    fclose(fimg);
    fclose(fcript);
    fclose(f);


}
void testul_chi_patrat(char *fisiersursa) {
    unsigned int* frr, *frg, *frb;
    double f, xr = 0, xg = 0, xb = 0, c;
    unsigned int h, w, i;

    FILE * fin = fopen(fisiersursa, "rb");

    if(fin==NULL) {

            printf("Fisierul nu a putut fi deschis\n");
            return;
    }
    fseek(fin, 18, SEEK_SET);
    fread(&w, sizeof(unsigned int), 1, fin);
    fread(&h, sizeof(unsigned int), 1, fin);
    pixel*img;
    img = (pixel*)malloc(h * w * sizeof(pixel));
    img = liniarizare(fisiersursa, &h, &w);

    frr = (unsigned int*)calloc(256, sizeof(unsigned int));
    frg = (unsigned int*)calloc(256, sizeof(unsigned int));
    frb = (unsigned int*)calloc(256, sizeof(unsigned int));

    for (i = 0; i < h * w; i++) {
            frr[img[i].r] = frr[img[i].r] + 1;
            frg[img[i].g] = frg[img[i].g] + 1;
            frb[img[i].b] = frb[img[i].b] + 1;
    }

    f = (h * w)/256;

    for (i = 0; i < 256; i++) {

       xr = xr + ((frr[i] - f) * (frr[i] - f)) / f;
       xg = xg + ((frg[i] - f) * (frg[i] - f)) / f;
       xb = xb + ((frb[i] - f) *(frb[i] - f)) / f;
    }

    printf("%.2f %.2f %.2f", xr, xg, xb);

}

void grayscale_image(char* nume_fisier_sursa, char* nume_fisier_destinatie) {
   FILE *fin, *fout;
   unsigned int dim_img, latime_img, inaltime_img;
   unsigned char pRGB[3], header[54], aux;

   printf("nume_fisier_sursa = %s \n", nume_fisier_sursa);

   fin = fopen(nume_fisier_sursa, "rb");
   if (fin == NULL) {
   		printf("nu am gasit imaginea sursa din care citesc");
   		return;
   	}

   fout = fopen(nume_fisier_destinatie, "wb+");

   fseek(fin, 2, SEEK_SET);
   fread(&dim_img, sizeof(unsigned int), 1, fin);
   printf("Dimensiunea imaginii in octeti: %u\n", dim_img);

   fseek(fin, 18, SEEK_SET);
   fread(&latime_img, sizeof(unsigned int), 1, fin);
   fread(&inaltime_img, sizeof(unsigned int), 1, fin);
   printf("Dimensiunea imaginii in pixeli (latime x inaltime): %u x %u\n", latime_img, inaltime_img);

   //copiaza octet cu octet imaginea initiala in cea noua
	fseek(fin, 0, SEEK_SET);
	unsigned char c;

	while (fread(&c, 1, 1, fin) == 1) {
		fwrite(&c,1,1,fout);
		fflush(fout);
	}

	fclose(fin);

	//calculam padding-ul pentru o linie
	int padding;
    if (latime_img % 4 != 0) {

            padding = 4 - (3 * latime_img) % 4;
    }
    else {
            padding = 0;
    }

    printf("padding = %d \n", padding);

	fseek(fout, 54, SEEK_SET);
	int i, j;
	for (i = 0; i < inaltime_img; i++) {
		for(j = 0; j < latime_img; j++) {
			//citesc culorile pixelului
			fread(pRGB, 3, 1, fout);
			//fac conversia in pixel gri
			aux = 0.299*pRGB[2] + 0.587*pRGB[1] + 0.114*pRGB[0];
			pRGB[0] = pRGB[1] = pRGB[2] = aux;
        	fseek(fout, -3, SEEK_CUR);
        	fwrite(pRGB, 3, 1, fout);
        	fflush(fout);
		}
		fseek(fout, padding, SEEK_CUR);
	}
	fclose(fout);
}
int **crearematrice(char*fisierimagine, unsigned int * h, unsigned int * w) {

    pixel **img;
    unsigned int i, j, inaltime_img, latime_img;
    FILE * fin = fopen(fisierimagine, "rb");

    if (fin == NULL) {
        printf("fisierul nu a putut fi deschis \n");
        return 0;
    }
    fseek(fin, 18, SEEK_SET);
    fread(&latime_img, sizeof(unsigned int), 1, fin);
    fread(&inaltime_img, sizeof(unsigned int), 1, fin);

    img = (pixel**)malloc(inaltime_img * (sizeof(pixel*)));

    for (i = 0; i < inaltime_img; i++) {

        img[i]=(pixel*)malloc(latime_img*(sizeof(pixel)));
    }

   fseek(fin, 54, SEEK_SET);
   int padding;
   if(latime_img % 4 != 0) {
        padding = 4 - (3 * latime_img) % 4;
   }
   else {

        padding = 0;
   }

   for(i = 0; i < inaltime_img; i++) {
        for(j = 0; j < latime_img; j++) {

           fread(&img[i][j].b,sizeof(unsigned char), 1, fin);
           fread(&img[i][j].g,sizeof(unsigned char), 1, fin);
           fread(&img[i][j].r,sizeof(unsigned char), 1, fin);

        }
        fseek(fin, padding, SEEK_CUR);
    }
    *h = inaltime_img;
    *w = latime_img;
    return img;


}

detectie* template_matching(pixel **img, unsigned int hi, unsigned int wi, char*fisiersablon, char*fisierdestinatiesab, double ps, unsigned int *nr) {

    grayscale_image(fisiersablon, fisierdestinatiesab);
    FILE *fsab = fopen(fisierdestinatiesab, "rb");

    if (fsab == NULL) {
            printf("Unul din fisiere nu poate fi deschis!\n");
            return;
    }

    unsigned int hs, ws, i, j;
    pixel **sablon;
    detectie *vect;
    sablon = crearematrice(fisierdestinatiesab, &hs, &ws);

    unsigned int cnt = 0, l, k, c = 0;
    double n;
    n = hs * ws;
    double smed = 0, fmed;

    for (i = 0; i < hs; i++) {
        for(j = 0; j < ws; j++) {
            smed=smed+(double)sablon[i][j].r;
        }
    }
    smed=smed/n;

    double sigmas = 0, sigmaf = 0, corr;

    for (i = 0; i < hs; i++) {
        for(j=0;j<ws;j++) {
                sigmas=sigmas+((double)sablon[i][j].r-smed)*((double)sablon[i][j].r-smed);
        }
    }
    sigmas = sigmas / (n - 1);
    sigmas = sqrt(sigmas);

    for (i = 0; i < hi; i++){
      for (j = 0; j < wi; j++) {
            if ((i + hs < hi) && (j + ws < wi)) {
                fmed = 0;
                sigmaf = 0;
                corr = 0;
                for (k = 0;k < hs; k++) {
                    for (l = 0; l < ws; l++){

                            fmed = fmed + img[i + k][j + l].r;
                    }
                }
                fmed = fmed / n;
                for (k = 0; k < hs; k++) {
                    for(l = 0; l < ws; l++) {
                            sigmaf = sigmaf + (img[i + k][j + l].r - fmed) * (img[i + k][j +l ].r - fmed);
                    }
                }
                sigmaf=sqrt(sigmaf/(n-1));

                for (k = 0; k < hs; k++) {
                    for (l = 0; l < ws; l++) {
                            corr = corr + (1 / (sigmaf * sigmas)) * (img[i + k][j + l].r - fmed) * (sablon[k][l].r - smed);
                    }
                }
                corr = corr / n;

                if (corr >= ps) {
                        if (cnt == 0) {
                            vect = (detectie*)malloc(sizeof(detectie));
                            vect[cnt].linie = i;
                            vect[cnt].coloana = j;
                            vect[cnt].corelatie = corr;
                            cnt = cnt + 1;
                        }
                        else {
                            vect = realloc(vect, (cnt + 1) * sizeof(detectie));
                            vect[cnt].linie = i;
                            vect[cnt].coloana = j;
                            vect[cnt].corelatie = corr;
                            cnt = cnt + 1;
                        }

             }


         }
      }
    }
     *nr = cnt;
     return vect;

}



void salveaza_imagine(pixel**img, unsigned int latime_img, unsigned int inaltime_img, char*nume_fisier_sursa, char* nume_fisier_destinatie) {

    unsigned int i;
    FILE *fin = fopen(nume_fisier_sursa, "rb");
    FILE *fout = fopen(nume_fisier_destinatie, "wb");

    if((fin == NULL) || (fout == NULL)) {
            printf("Fisierul in care dorim sa salvam extern nu poate fi deschis!\n");
            return;
    }
    unsigned char c;
    while (fread(&c, 1, 1,fin) == 1) {

		fwrite(&c, 1, 1,fout);
		fflush(fout);

	}
	int padding;
    if (latime_img % 4 != 0) {
         padding = 4 - (3 * latime_img) % 4;
    }
    else {
         padding = 0;
    }
    unsigned int j;
	fseek(fout, 54, SEEK_SET);

    for (i = 0; i < inaltime_img; i++){
            for (j = 0; j < latime_img; j++) {

                fwrite(&img[i][j].b, sizeof(unsigned char), 1, fout);
                fflush(fout);
                fwrite(&img[i][j].g, sizeof(unsigned char), 1, fout);
                fflush(fout);
                fwrite(&img[i][j].r, sizeof(unsigned char), 1, fout);
                fflush(fout);
            }
            fseek(fout, padding, SEEK_CUR);
    }
    fclose(fin);
    fclose(fout);
}


pixel** contur(pixel**img, unsigned int hi, unsigned int wi, detectie fereastra, pixel culoare, unsigned int hs, unsigned int ws) {
    unsigned int i,j;

    for (j = 0; j < ws; j++) {

         img[fereastra.linie][fereastra.coloana+j].r = fereastra.culoare.r;
         img[fereastra.linie][fereastra.coloana+j].g = fereastra.culoare.g;
         img[fereastra.linie][fereastra.coloana+j].b = fereastra.culoare.b;
         img[fereastra.linie + hs - 1][j + fereastra.coloana].r = fereastra.culoare.r;
         img[fereastra.linie + hs - 1][j + fereastra.coloana].g = fereastra.culoare.g;
         img[fereastra.linie + hs - 1][j + fereastra.coloana].b = fereastra.culoare.b;
    }

    for (i = 0; i < hs; i++){

         img[fereastra.linie + i][fereastra.coloana].r = fereastra.culoare.r;
         img[fereastra.linie + i][fereastra.coloana].g=fereastra.culoare.g;
         img[fereastra.linie + i][fereastra.coloana].b=fereastra.culoare.b;
         img[fereastra.linie + i][fereastra.coloana + ws - 1].r = fereastra.culoare.r;
         img[fereastra.linie + i][fereastra.coloana + ws - 1].g = fereastra.culoare.g;
         img[fereastra.linie + i][fereastra.coloana + ws - 1].b = fereastra.culoare.b;
    }

    return img;


}

int cmp(const void*a, const void*b) {
    double l, r;
    l = ((detectie*)a) -> corelatie;
    r = ((detectie*)b) -> corelatie;
    if (r - l > 0) {
            return 1;
    }
    else {
            return -1;
    }
}


detectie* ordonare(detectie *d, unsigned int k) {
    qsort(d, k, sizeof(detectie), cmp);
    return d;

}


int suprapunere(detectie di, detectie dj, unsigned int h, unsigned int w) {

    double ariedi, ariedj, intersectie;
    unsigned int hi = 0, wi = 0;
    ariedi = w * h;
    ariedj = w * h;

    if ((di.linie <= dj.linie) && (dj.linie <= (di.linie + h - 1))) {
            hi = h - (dj.linie - di.linie);
            if ((di.coloana <= dj.coloana) && (dj.coloana <= (di.coloana + w - 1))) {
                    wi=w-(dj.coloana-di.coloana);
            }
            else {
                    if((di.coloana>=dj.coloana)&&(dj.coloana+w-1>=di.coloana)){
                            wi=w-(di.coloana-dj.coloana);
                    }
                    else {
                            return 0;
                    }
            }
    }



    else {
            if ((di.linie >= dj.linie) && (di.linie <= (dj.linie + h - 1))) {

                hi = h - (di.linie - dj.linie);
                if ((di.coloana <= dj.coloana) && (dj.coloana <= (di.coloana + w - 1))) {
                        wi = w - (dj.coloana - di.coloana);
                }
                else {
                        if ((di.coloana >= dj.coloana) && (dj.coloana + w - 1 >= di.coloana)) {

                            wi=w-(di.coloana-dj.coloana);
                        }
                        else {
                                    return 0;
                             }
                }
            }
            else {
                    return 0;
            }

    }

    intersectie = (hi * wi) / (ariedi + ariedj - (hi * wi));

    if (intersectie > 0.2) {
            return 1;
    }
    else {
            return 0;
    }
}

detectie *elimin_non_maxime(detectie *d, unsigned int k, unsigned int h, unsigned int w, unsigned int *nr) {
    unsigned int *vizitat;
    unsigned int cnt = 0, i, j, l;
    detectie *vector;
    vizitat = (unsigned int*)calloc(k, sizeof(unsigned int));

    for (i = 0; i < k; i++) {
        for ( j = i + 1; j < k; j++) {
                if (vizitat[j] == 0) {
                    vizitat[j] = suprapunere(d[i], d[j], h, w);
                }
        }
    }

    for (i = 0; i < k; i++) {

            if(vizitat[i] == 0) {
                    cnt = cnt + 1;
            }
    }
    vector = (detectie*)malloc(cnt * sizeof(detectie));
    l = 0;

    for (i = 0; i < k; i++) {
            if (vizitat[i] == 0) {
                vector[l] = d[i];
                l = l + 1;
            }
    }
    *nr=cnt;

    return vector;


}


int main() {

   char img[35], img_criptata[35], secret_key[35];

   printf("Numele fisierului care contine imaginea care va fi criptata: ");
   fgets(img, 35, stdin);
   img[strlen(img) - 1] = '\0';

   printf("Numele fisierului unde va fi salvata criptarea : ");
   fgets(img_criptata, 35, stdin);
   img_criptata[strlen(img_criptata) - 1] = '\0';

   printf("Numele fisiserului din care se citeste cheia : ");
   fgets(secret_key, 35, stdin);
   secret_key[strlen(secret_key) - 1] = '\0';

   criptare(img, img_criptata, secret_key);

   char img_decriptata[35], img_criptata2[35], cheiasecreta[35];

   printf("Numele imaginii in care se va salva decriptarea : ");
   fgets(img_decriptata, 35, stdin);
   img_decriptata[strlen(img_decriptata) - 1] = '\0';

   printf("Numele imaginii care se decripteaza : ");
   fgets(img_criptata2, 35, stdin);
   img_criptata2[strlen(img_criptata2) - 1] = '\0';

   printf("Numele fisierului din care se citeste cheia secreta : ");
   fgets(cheiasecreta, 35, stdin);
   cheiasecreta[strlen(cheiasecreta) - 1] = '\0';

   decriptare(img_decriptata,img_criptata2,cheiasecreta);
   printf("\n");

   printf("Testul chi_patrat pt imaginea initiala: \n");
   testul_chi_patrat(img);
   printf("\n");

   printf("testul chi_patrat pt imaginea criptata: \n");
   testul_chi_patrat(img_criptata);
   printf("\n");

   double ps=0.5;
   unsigned int nrsab;

   printf("Se citeste numarul de sabloane : ");
   scanf("%u", &nrsab);
   fgetc(stdin);

   detectie **v;
   pixel *culoare;

   culoare = (pixel*)malloc(nrsab * sizeof(pixel));
   culoare[0].r = 255;
   culoare[0].g = 0;
   culoare[0].b = 0;
   culoare[1].r = 255;
   culoare[1].g = 255;
   culoare[1].b = 0;
   culoare[2].r = 0;
   culoare[2].g = 255;
   culoare[2].b = 0;
   culoare[3].r = 0;
   culoare[3].g = 255;
   culoare[3].b = 255;
   culoare[4].r = 255;
   culoare[4].g = 0;
   culoare[4].b = 255;
   culoare[5].r = 0;
   culoare[5].g = 0;
   culoare[5].b = 255;
   culoare[6].r = 192;
   culoare[6].g = 192;
   culoare[6].b = 192;
   culoare[7].r = 255;
   culoare[7].g = 140;
   culoare[7].b = 0;
   culoare[8].r = 128;
   culoare[8].g = 0;
   culoare[8].b = 128;
   culoare[9].r = 128;
   culoare[9].g = 0;
   culoare[9].b = 0;

   unsigned int *nr, wi, hi, hs, ws;
   pixel **imag, **sab;
   unsigned int j, i;
   char imagine[35];

   printf("Numele imaginii pe care se va rula template matching : ");
   fgets(imagine, 35, stdin);
   imagine[strlen(imagine) - 1] = '\0';
   printf("\n");

   char sablon[30][35];

   for (i = 0; i < nrsab;i++) {
        printf("Se citesteste numele sablonului cu nr %u: ", i);
        fgets(sablon[i], 35, stdin);
        sablon[i][strlen(sablon[i]) - 1] = '\0';
        printf("\n");

    }

   grayscale_image(imagine, "testgary.bmp");
   imag = crearematrice("testgary.bmp", &hi, &wi);

   v = (detectie**)malloc(nrsab * sizeof(detectie*));
   nr = (unsigned int)malloc(nrsab * sizeof(unsigned int));

   for (i = 0; i < nrsab; i++) {
        v[i]=template_matching(imag,hi,wi,sablon[i],"cifra0gray.bmp",ps,&nr[i]);
   }

    detectie *d;
    unsigned int k=0;

    for (i = 0; i < nrsab; i++) {
            k = k + nr[i];
    }

    d = (detectie*)malloc(k * sizeof(detectie));
    unsigned int p = 0;

    for (i = 0; i < nrsab; i++) {
        for (j = 0; j < nr[i]; j++) {
                d[p].linie = v[i][j].linie;
                d[p].coloana = v[i][j].coloana;
                d[p].corelatie = v[i][j].corelatie;
                d[p].culoare.r = culoare[i].r;
                d[p].culoare.g = culoare[i].g;
                d[p].culoare.b = culoare[i].b;
                p=p+1;
        }
    }

    detectie *df;

    df = ordonare(d, k);

    detectie *fin;
    pixel**imagnormala;

    imagnormala = crearematrice(imagine, &hi, &wi);

    unsigned int cnt;

    sab = crearematrice(sablon[0], &hs, &ws);

    fin = elimin_non_maxime(df, k, hs, ws, &cnt);

    for (i = 0; i < cnt; i++){
            imagnormala = contur(imagnormala , hi , wi , fin[i] , culoare[i] , hs,ws);
    }

    char imaginedupacontur[35];
    printf("\n");

    printf("Se citeste numele fisierului in care se va salva imaginea dupa ce vor fi desenate contururile : ");
    fgets(imaginedupacontur, 35, stdin);
    imaginedupacontur[strlen(imaginedupacontur) - 1] = '\0';

    salveaza_imagine(imagnormala, wi, hi, imagine, imaginedupacontur);
    return 0;
}
