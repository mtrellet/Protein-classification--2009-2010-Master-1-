# -*- coding: cp1252 -*-
from Bio import ExPASy
from Bio import SeqIO
import webbrowser
from Tkinter import *
from ScrolledText import ScrolledText
import tkFileDialog
import tkFont
import tkMessageBox
import time as origtime

                #################################################################
                #                                                               #
                #       R�cup�ration du chemin d'acc�s au fichier m8,           #
                #       des e-value et pourcentage d'identit� saisis            #
                #       par l'utilisateur                                       #
                #                                                               #
                #       Auteurs : Elise Larsonneur, Mika�l Trellet              #
                #                                                               #
                #################################################################

def verif_value():
        #Recup�ration de la cha�ne de caract�res correspondant au nom de fichier ins�r�
       #dans la zone de texte 'entrPar' d'insertion du nom de fichier de donn�es blastall apr�s clic sur le bouton 'Parcourir'
        RecupFich=entrPar.get()
        RecupEval=entr2a.get()
        RecupPid=entr2b.get()
        etat1=""
        etat2=""
        
        if(RecupFich==''):
                tkMessageBox.showinfo("Erreur","Nom de fichier non sp�cifi�")
                etat="non"
        else:
                
                try:
                        if ((type(eval(RecupEval)) == float) and (eval(RecupEval)>0) and (eval(RecupEval)<1)):
                                etat1="ok"
                
                        else:
                                raise TypeError
                        
                except (SyntaxError,TypeError, NameError,ValueError):
                        tkMessageBox.showinfo("Erreur","Entrer une e-value de type r�el positif")
                        etat1="non"
                
                try:
                        if ((type(eval(RecupPid ))== int) and (eval(RecupPid)>0) and (eval(RecupPid)<=100)):
                                etat2="ok"
                
                        else:
                                raise TypeError
                       
                except (SyntaxError,TypeError, NameError,ValueError):
                        tkMessageBox.showinfo("Erreur","Entrer un pourcentage d'identit� de type r�el positif")
                        etat2="non"

        print etat1, etat2

        if((etat1=="ok") and (etat2=="ok")):
                ObFichier = open(RecupFich,'r')
                lancer_test(ObFichier, RecupEval, RecupPid)

def lancer_test(ObFichier, RecupEval, RecupPid):

          ########################################
          #                                      #
          # R�cup�ration des donn�es BLASTall    #
          #                                      #
          #  Auteur : Amyra Aliouat              #
          #                                      #  
          #       Optimisation : Mika�l Trellet  #
          #                                      #
          ########################################

        ##### Changer de r�pertoire courant #########
        #from os import chdir
        #chdir("C:\Users\selfharm\Desktop\PutaindPython\Projet python")
                
        debut=origtime.time()

        ####### Ouverture du fichier, initialisation d'une liste pour r�cup�rer les donn�es, du dictionnaire principal et des variables n�cessaires aux traitements des prochaines �tapes. #######
        #ObFichier = open('Paramecie_test.blastall.m8','r')
        #ObFichier = open('Ptetraurelia_peptides_cur.blastall.m8','r')
        
        liste = []
        Dico={}
        length_tot = 0
        length_average = 0.0
        nb_genes = 0
        i=0
        precedent=""
        percentid=eval(RecupPid)
        evalue=eval(RecupEval)
        ##### lecture et stockage de toutes les lignes du fichier ########

        for line in ObFichier.xreadlines():
                i=i+1
                if (i%100000==0):
                        temps=origtime.time()-debut
                        
                liste = []
                liste = line.split('\t')
                if(len(liste)==12):
                        #liste[0]:query id, liste[1]:subject id, liste[2]:%identity, liste[3]:alignment length, liste[6]:query start, liste[7]:query end, 
#liste[8]:subject start, liste[9]:subject end, liste[10]:expect value.
                        liste[10]= float(liste[10])
                        liste[2]=float(liste[2])

                        if ((precedent==liste[0]) and (liste[0]!=liste[1]) and (liste[10]<evalue) and (liste[2]>percentid)):   
                                sousDicoQuery={}
                                Dico[liste[0]][liste[1]]=0
                                sousDicoQuery[liste[1]]=liste[2],liste[3],liste[6],liste[7],liste[8],liste[9],liste[10]
                                Dico[liste[0]][liste[1]]=sousDicoQuery.items()[0][1]
       
                        elif ((precedent!=liste[0]) and (liste[0]!=liste[1]) and (liste[10]<evalue) and (liste[2]>percentid)):
                                Dico[liste[0]]=0
                                sousDicoQuery={}
                                sousDicoQuery[liste[1]]=liste[2],liste[3],liste[6],liste[7],liste[8],liste[9],liste[10]
                                Dico[liste[0]]=sousDicoQuery#.items()
                                precedent=liste[0]

                        if(liste[0]==liste[1]):
                                nb_genes = int(nb_genes+1)
                                length_tot = length_tot + int(liste[3])
                        
        
        ObFichier.close()
        keys=Dico.keys()
        print "nb de genes :",nb_genes
        print "long tot :", length_tot
        length_average = float((length_tot)/(nb_genes))
        print "Longueur moyenne des genes :", length_average
        print "nb queries gardes =",len(Dico)
        temps=origtime.time()-debut
        print ("temps : "+str(temps))
        print ("nbre de lignes ds fichier : "+str(i))


        Dico_save={}
        Dico_save.update(Dico);

        print "Dico_save : ", len(Dico_save)


                #################################################################
                #                                                               #
                #                       Clustering                              #
                #                                                               #
                #                  Auteur : Benjamin Zerath,                    #
                #               Optimisation : Mika�l Trellet                   #
                #                                                               #
                #################################################################


        def clustering(Dico):
                print len(Dico)
                j=0
                nbClef=0
       
                for key in keys:#parcourt le grand dico ex : 160001
                        if (j%10000==0):
                                #print "10000 : ",j
                                temps=origtime.time()-debut
                                #print temps
                        j=j+1
                        #print "key1 :",key
                        if (Dico.has_key(key)):
                                sousKeysDico=Dico[key].keys()
                                #print "avant",Dico[key]
                #        print Dico[souskey]
                                for souskey in sousKeysDico: #parcourt le sous-dico ex : 366001
                                        if (Dico.has_key(souskey)): #v�rifie que les subjects sont aussi l� en tant que query qq part
                                                souskeysSubject=Dico[souskey].keys()#prend les cl�s des subject du subject en cours
                                                for souskeysubject in souskeysSubject:
                                                        if (souskeysubject not in sousKeysDico):
                                                                Dico[key][souskeysubject]=Dico[souskey][souskeysubject]
                                                                
                                                        nbClef=len(Dico[souskey].keys())
                                                               
                                if (souskey!=key) and (Dico.has_key(souskey)):
                                        #print "sous key : ",souskey
                                        del Dico[souskey]
                                #print "apr�s",Dico[key]
                print "\n \n \n"
                print "longueur dico", len(Dico)
                

                temps2=origtime.time()-debut
                print temps2


           


        ######################################################################################
        #                                                   				     #
        #                 Attribution de fonctions aux prot�ines :			     #
        #										     #
        #Associe � chaque prot�ine sa fiche swiss Prot et ses identifiants de Gene Ontology  #
        #    						    				     #
        #                 Auteur: Marie       le 04.01.2010              		     #
        #                                                   				     #
        #          Repertoire: /ProjectPythonFonction/Python_finale                          #
        # 										     #
        ######################################################################################
 

        ##### I Creation du dictionnaire Dico {id_GO : fonction precise , fonction generale} ######

        ####### Stockage de toute les lignes du fichier 'process_bio' dans la variable lecture #######
        objF1 = open('process_bio.txt','r')
        lecture = objF1.readlines()


        ####### Recuperation des lignes commencant par 'id' et 'name' et  stockage dans les variable liste1,liste2 ######
        liste1=[]
        liste2=[]
        ide="id"
        fct="name"
        fctg="biological process"

        for line in lecture:            ## pour chaque ligne du fichier lecture
          if line.startswith(ide):      ## si une ligne commence par ide ("id")
            iden=line
            iden=iden.strip()           ## supprime le saut de ligne \n a la fin de iden
            iden=iden.replace('id: ','')## remplace id: par un espace
            liste1.append(iden)         ## met tous les iden dans liste1
   

          if line.startswith(fct):          ## si une ligne commence par fct ("name")
            fonct=line
            fonct=fonct.strip()             ## supprime le saut de ligne \n a la fin de fonct
            fonct=fonct.replace('name: ','')## remplace name: par un espace
            liste2.append(fonct)            ## met tous les fonct dans liste2
        objF1.close()                       ## Fermeture du fichier objF1

        ###### Stockage de liste1,liste2 dans Dico ######
        Dico1={}
        for i in range(len(liste1)):        ## pour chaque element i de la liste1 
          Dico1[liste1[i]]=liste2[i],fctg    ## on cr�� un Dico1 avec la liste1 comme cl� et on associe 2 valeurs au Dico1 : liste2 et fctg [fonction generale]



        ####### Stockage de toute les lignes du fichier compo_cell dans la variable lecture #######
        objF2 = open('compo_cell.txt','r')
        lecture = objF2.readlines()

        ####### Recuperation des lignes commencant par 'id' et 'name' et stockage dans liste1,liste2 ######
        liste1=[]
        liste2=[]
        ide="id"
        fct="name"
        fctg="cellular component"

        for line in lecture:              ## pour chaque ligne du fichier lecture
          if line.startswith(ide):        ## si une ligne commence par ide ("id")
            iden=line
            iden=iden.strip()             ## supprime le saut de ligne \n a la fin de iden
            iden=iden.replace('id: ','')  ## remplace id: par un espace
            liste1.append(iden)           ## met tous les iden dans liste1
   

          if line.startswith(fct):          ## si une ligne commence par fct ("name")
            fonct=line
            fonct=fonct.strip()             ## supprime le saut de ligne \n a la fin de fonct
            fonct=fonct.replace('name: ','')## remplace name: par un espace
            liste2.append(fonct)            ## met tous les fonct dans liste2
        objF2.close()                       ## Fermeture du fichier objF2


        ###### Stockage de liste1,liste2 dans Dico1 ######
        for i in range(len(liste1)):        ## pour chaque element i de la liste1 
          Dico1[liste1[i]]=liste2[i],fctg    ## on cr�� un Dico1 avec la liste1 comme cl� et on associe 2 valeurs au Dico1 : liste2 et fctg [fonction generale]
  


        ####### Stockage de toute les lignes du fichier fonc_mol dans la variable lecture #######
        objF3 = open('fon_mol.txt','r')
        lecture = objF3.readlines()

        ####### Recuperation des lignes commencant par 'id' et 'name' et stockage dans des listes ######
        liste1=[]
        liste2=[]
        ide="id"
        fct="name"
        fctg="molecular fonction"

        for line in lecture:              ## pour chaque ligne du fichier lecture
          if line.startswith(ide):        ## si une ligne commence par ide ("id")
            iden=line
            iden=iden.strip()             ## supprime le saut de ligne \n a la fin de iden
            iden=iden.replace('id: ','')  ## remplace id: par un espace
            liste1.append(iden)           ## met tous les iden dans liste1
   

          if line.startswith(fct):          ## si une ligne commence par fct ("name")
            fonct=line                      
            fonct=fonct.strip()             ## supprime le saut de ligne \n a la fin de fonct
            fonct=fonct.replace('name: ','')## remplace name: par un espace
            liste2.append(fonct)            ## met tous les fonct dans liste2
        objF3.close()                       ## Fermeture du fichier objF3

        ###### Stockage de liste1,liste2 dans Dico1 ######
        for i in range(len(liste1)):        ## pour chaque element i de la liste1 
          Dico1[liste1[i]]=liste2[i],fctg    ## on cr�� un Dico1 avec la liste1 comme cl� et on associe 2 valeurs au Dico1 : liste2 et fctg [fonction generale]

       
        ########################################################################################


        ##### II Creation du dictionnaire SW{idPDB:idSW} (idPDB=id_paramecium ; idSW=lien vers fiche SW) ######

        ####### Stockage du numero d'accession du gene et sa fiche swiss prot dans un dictionnaire a partir du fichier protAC #######
        ObFichier = open('protAC.txt','r')
        ligne = ObFichier.readline()

        ###### stockage du numero d'accession et la fiche swiss prot sous forme de liste  ####
        SW={}
        liste1=[]
        liste2=[]
        liste =[]
        liste = ligne.split('\t')   ## Decoupage de chaque ligne du fichier en liste apres une tabulation

        ## Recuperation de la 1ere ligne du fichier
        idPDB = liste[0]            ## Stocke la valeur de liste[0] dans idPDB
        liste1.append(idPDB)        ## met tous les idPDB dans liste1
        idSW = liste[3]             ## Stocke valeur de liste[3] dans idSW
        idSW=idSW.strip()           ## supprime le saut de ligne \n a la fin de idSW
        liste2.append(idSW)         ## met tous les idSW dans liste2
        SW[idPDB]=idSW              ## Dico1 SW     

        #### stockage du numero d'accession et la fiche swiss prot dans un dictionnaire SW{idPDB:idSW}  ######
        for ligne in ObFichier:     ## pour chaque ligne du fichier ObFichier
          liste = ligne.split('\t') ## Decoupage de chaque ligne du fichier en liste apres une tabulation
          if (liste[0]!=''):
            idPDB = liste[0]          ## Stocke valeur de liste[0] dans idPDB
            liste1.append(idPDB)      ## met tous les idPDB dans liste1
            idSW = liste[3]           ## Stocke valeur de liste[3] dans idSW
            idSW=idSW.strip()         ## supprime le saut de ligne \n a la fin de idSW
            liste2.append(idSW)       ## met tous les idSW dans liste2
            SW[idPDB]=idSW            ## Dico1 SW
        ObFichier.close()           ## Fermeture du fichier ObFichier

      
        ##### III Creation du dictionnaire Pinf{idPDB:fonction,deb,fin,scaffold} ######

        ####### Stockage du numero d'accession du gene et sa fonction dans un dictionnaire a partir du fichier protein_info #######
        ObFichier = open('protein_info.txt','r')
        ligne = ObFichier.readline()

        ###### Stockage du numero d'accession et sa fonction sous forme de liste ####
        Pinf={}
        liste1=[]
        liste2=[]
        taille_tot=0 ## Creation d'une variable pour stocker la taille de chaque gene, en vue de faire une moyenne

        #### Creation du dictionnaire Pinf{idPDB:fonction,deb,fin,scaffold}####
        for ligne in ObFichier:                         ## pour chaque ligne du fichier ObFichier
          liste = []
          liste = ligne.split('\t')                     ## Decoupage de chaque ligne du fichier en liste apres une tabulation
          if(liste[0]!=''):                             ## pour chaque num�ro d'accession (Si la 1er colonne n'est pas vide)
            idPDB = liste[0]                            ## Stocke les valeurs de liste[0] dans idPDB
            liste1.append(idPDB)                        ## met tous les idPDB dans liste1
            deb = liste[3]                              ## Stocke les valeurs de liste[3] dans deb
            fin = liste[4]                              ## Stocke les valeurs de liste[4] dans fin
            scaffold = liste[2]                         ## Stocke les valeurs de liste[2] dans scaffold

        #### lorsqu'aucune fonction n'est pas associ� aun gee, on attrubue :'no_function' au gee 
            if(liste[6] != ''):                         ## Si la 5eme colonne n'est pas vide (existence d'une fonction)
              finfo = liste[6]                          ## Stocke valeur de liste[6] dans finfo
            else:                                       ## il n'y a pas de fonction associ� a un gene (liste[6] est vide)
              finfo = 'no_function'                     ## on ecrit 'no_function' pour les gene qui n'ont pas de fonction
            finfo=finfo.strip()                         ## supprime le saut de ligne \n a la fin de finfo
            liste2.append(finfo)                        ## met tous les finfo dans liste2
            Pinf[idPDB]=finfo, deb, fin, scaffold       ## Dico1 Pinf
            taille_tot+=int(fin) - int(deb)
        ObFichier.close()                               ## Fermeture du fichier ObFichier  

        
           ################################################################################

                #################################################################
                #                                                               #
                #       Cr�ation Dictionnaire DicoGO[idPDB]:idGO                #
                #                                                               #
                #                  Auteur : Elise Larsonneur                    #
                #                                                               #
                #################################################################

        Obfichier=open('Liste_idPDB_idGO.txt','r')
        DicoGO={}
        for line in Obfichier.xreadlines():
                liste = []
                liste = line.split('\t')
                liste1 =liste[2].replace('[','')
                liste2 =liste1.replace(']','')
                liste3 =liste2.replace("'",'')
                liste4 =liste3.replace('\n','')
                liste5 =liste4.replace(' ','')
                liste6 =liste5.split(",")
                DicoGO[liste[1]]=liste6
        
        Obfichier.close()
        print "DicoGO : ",len(DicoGO)

        ################################################################################


                #################################################################
                #                                                               #
                #   Distinction Duplications Tandem - Duplications Segmentales  #
                #                                                               #
                #                  Auteur : Mika�l Trellet                      #
                #                                                               #
                #################################################################
                
        taille_moyenne=(taille_tot/len(Pinf))

        print "taille moyenne =",taille_moyenne

        print "Pinf ",len(Pinf)
        print "SW ",len(SW)
        print "Dico1 ",len(Dico1)
        print "Dico_save : ",len(Dico_save)


        Tandem={}
        Segmental={}
        i=0
        for key in Dico_save:   ##pour chaque prot�ine cl� du dictionnaire Dico-save
                if (Pinf.has_key(key)):        ##si cette prot�ine est contenue dans le dictionnaire Pinf
                        sousDicoquery1={}
                        sousDicoquery2={}
                        souskeys=Dico_save[key].keys();    ##le nouveau dictionnaire souskeys contient les prot�ines sous-cl�s de Dico_save
                
                        for souskey in souskeys:
                                #####si le couple de prot�ines cl�-sous cl� est sur le m�me scaffold (m�me chromosome) 
                                if (Pinf.has_key(souskey)) and (Pinf[key][3]==Pinf[souskey][3]):
                                        ##si la position de d�but de la prot�ine cl� est sup�rieure � la position de fin de la prot�ine sous-cl�
                                        ##et la distance entre la position de d�but de la prot�ine cl� et la position de fin de la prot�ine sous-cl�
                                        ##est inf�rieure � la taille moyenne d'une prot�ine
                                        ##ou si la position de d�but de la prot�ine cl� est sup�rieure � la position de fin de la prot�ine sous-cl�
                                        ##et la distance entre la position de d�but de la prot�ine sous-cl� et la position de fin de la prot�ine cl�
                                        ##est inf�rieure � la taille moyenne d'une prot�ine
                                        if((int(Pinf[key][1])>int(Pinf[souskey][2])) and ((int(Pinf[key][1])-int(Pinf[souskey][2]))<taille_moyenne))or ((int(Pinf[souskey][1])>int(Pinf[key][2])) and (int(Pinf[souskey][1]) - int(Pinf[key][2])<taille_moyenne)):
                                                ##alors le dictionnaire suivant contient ce couple de prot�ines  (duplication tandem)
                                                sousDicoquery1[souskey]=Dico[key][souskey]
                                        else:
                                                ##sinon le 2nd dictionnaire contient ce couple de prot�ines (duplication segmentale)
                                                sousDicoquery2[souskey]=Dico[key][souskey]
                                                
                                ##si le couple de prot�ines cl�-sous cl� n'est pas sur le m�me scaffold 
                                else:
                                        if (Pinf.has_key(souskey)):
                                                ##le 2nd dictionnaire contient ce couple de prot�ines (duplication segmentale)
                                                sousDicoquery2[souskey]=Dico[key][souskey]
                                                
                
                        if(len(sousDicoquery1)!=0):
                                Tandem[key]=sousDicoquery1
                        if(len(sousDicoquery2)!=0):
                                Segmental[key]=sousDicoquery2
                                
                ## si le dictionnaire Pinf ne contient pas une cl� de Dico_save et le couple de prot�ines non obligatoirement sur le m�me scaffold :                
                else:
                        i+=1
                        Segmental[key]=Dico[key]   ##le couple cl�-souscl� est consid�r�e comme duplication segmentale
                                

        print " i : ",i        
        print "Tandem : ", len(Tandem)
        print "Segmental : ", len(Segmental)

        clustering(Tandem)
        clustering(Segmental)
        clustering(Dico_save)
        

        #################################################################################
        #                                                                               #
        #         Interface, affichage des r�sultats, "Lien" vers Fiche Uniprot         #
        #                                                                               #
        #                       Auteur : Elise Larsonneur                               #
        #                                                                               #
        #################################################################################

        
        #############Affichage des r�sultats dans une fen�tre############
    


        #######Fonction permettant d'enregistrer dans un fichier la liste de duplications tandem######
        def enregistrerDupliT():
            filename = tkFileDialog.asksaveasfilename(filetypes = [('Fichiers texte', '*.txt'),('Fichiers doc','*.doc'),('Fichiers odt','*.odt')])
    
            if (not filename):
                tkMessageBox.showinfo("Erreur","Nom de fichier non sp�cifi�")
            else :
                fichEnreg=open (filename,'a')
                string=(txt.get(1.0,END))
                fichEnreg.write(string.encode('utf-8'))#encode() convertit un objet Unicode en objet string 8bit
                fichEnreg.close()
        #######Fonction permettant d'enregistrer dans un fichier la liste de duplications segmentales######

        def enregistrerDupliS():
            filename = tkFileDialog.asksaveasfilename(filetypes = [('Fichiers texte', '*.txt'),('Fichiers doc','*.doc')('Fichiers odt','*.odt')])
    
            if (not filename):
                  tkMessageBox.showinfo("Erreur","Nom de fichier non sp�cifi�")
            else :
                fichEnreg=open (filename,'a')
                string=(txt3.get(1.0,END))
                fichEnreg.write(string.encode('utf-8'))#encode() convertit un objet Unicode en objet string 8bit
                fichEnreg.close()


        ############Creation de fen�tres contenant une zone de texte associ�e � une barre de d�filement#########
        
        #########Cr�ation d'une premi�re fen�tre pour l'insertion d'une liste de duplications tandem##########
        fen2=Tk()
        fen2.title('Duplications Tandem')
        scroll = Scrollbar(fen2)
        txt = Text(fen2,width=150,height=30)
        txt.focus_set()
        scroll.pack(side=RIGHT, fill=Y)
        txt.pack(side=LEFT, fill=Y)
        scroll.config(command=txt.yview)
        txt.config(yscrollcommand=scroll.set)

        mainmenu = Menu(fen2)  ###Cr�ation d'un menu Fichier
        menuFichier = Menu(mainmenu)
        menuFichier.add_command(label="Enregistrer sous",command=enregistrerDupliT)  
        menuFichier.add_separator() ## Ajout d'une ligne s�paratrice
        menuFichier.add_command(label="Quitter", command=fen2.quit)
        mainmenu.add_cascade(label = "Fichier", menu = menuFichier)
        fen2.config(menu = mainmenu)


        ############Cr�ation d'une premi�re fen�tre pour l'insertion d'une liste de duplications segmentales######
        fen3=Tk()
        fen3.title('Duplications Segmentales')
        scroll3 = Scrollbar(fen3)
        txt3 = Text(fen3,width=200,height=30)
        txt3.focus_set()
        scroll3.pack(side=RIGHT, fill=Y)
        txt3.pack(side=LEFT, fill=Y)
        scroll3.config(command=txt3.yview)
        txt3.config(yscrollcommand=scroll3.set)
        
        mainmenu2 = Menu(fen3)  ###Cr�ation d'un menu Fichier
        menuFichier2 = Menu(mainmenu2)
        menuFichier2.add_command(label="Enregistrer sous",command=enregistrerDupliS)  
        menuFichier2.add_separator() ## Ajout d'une ligne s�paratrice
        menuFichier2.add_command(label="Quitter", command=fen3.quit)
        mainmenu2.add_cascade(label = "Fichier", menu = menuFichier2)
        fen3.config(menu = mainmenu2)


        #########affichage des info d'une proteine###################

        def affichageInfoGene(idPDB,txt):

          #recuperation des idGO de la proteine identifiee par idPDB sur le site Uniprot  
          if (Pinf.has_key(idPDB)):    
                debut=Pinf.get(idPDB)[1]  #position de d�but sur le chromosome
                fin=Pinf.get(idPDB)[2]    #position de fin sur le chromosome
                chrom=Pinf.get(idPDB)[3]   #scaffold : chromosome
                fctProtInf=Pinf.get(idPDB)[0]   #fonction tr�s pr�cise de la proteine (fct protinfo)
                idGO=DicoGO.get(idPDB)   #liste d'idGO si existence d'idGO dans Uniprot, "noGO" sinon
                
          else:
               debut='unknown'
               fin='unknown'
               chrom='unknown'
               fctProtInf='unknown'
               idGO='unknown'

           ##insertion des donnees dans la zone de texte de la fen�tre
               
          string=idPDB+" "+debut+" "+fin+" "+chrom#+"\n"
          txt.insert(END, string)
               
          if (idGO[0])!=('noGO'):     #si la proteine possede un ou plusieurs idGO :
                for i in range (0,len(idGO)):    #pour chaque idGO de la liste d'idGO de la proteine
                        if ((idGO=='unknown')):
                                string1="\n\t"+idGO+"//" +fctProtInf
                        elif (idGO!='unknown') and (not Dico1.get(idGO[i])):
                                string1="\n\t"+idGO[i]+"//" +fctProtInf
                        else:
                                string1="\n\t"+idGO[i]+" // "+Dico1.get(idGO[i])[1] +" // "+ Dico1.get(idGO[i])[0]+" // "+fctProtInf
                                        

                        txt.insert(END, string1)  #insertion de la cha�ne de caract�res dans la zone de texte de la fen�tre
                
          else:
                    #si la proteine ne possede pas d'idGO sur le site Uniprot :
                    string2= "\n\t"+idGO[0]+" // "+fctProtInf  #idGO, fct protinfo
                    txt.insert(END, string2)
            
          txt.insert(END, "\n")
        
           

        #tests
        #affichageInfoGene('GSPATP00037024001')  #Proteine "segmental"  ss IdGO
        #affichageInfoGene('GSPATP00012642001')  #Proteine "segmental" ac IdGO

        ###########################################


        #######Fonction permettant l'affichage d'une liste de duplications dans la fen�tre correspondante######

        def affichageCluster(Dico,fen) :  #prend en entr�e un dictionnaire et une fen�tre
        
                j=0
                print len(Dico)
                for key in Dico:
                    
                    if (Dico.has_key(key)):
                          j=j+1
                          sousKeysDico=Dico[key].keys()

                          txt.insert(END, ("Cluster "+str(j)+"\n"))
                          txt.insert(END, "Prot�ine r�f�rence: ")
                          affichageInfoGene(str(key),txt)
                          txt.insert(END,"\n")
                          #print Dico[key].keys()
                          #print Dico[key]
                          
                
                          for souskey in sousKeysDico:
                                  if (souskey!=key):
                                          affichageInfoGene(souskey,txt)
                                          txt.insert(END,"\n")
                          txt.insert(END,"\n\n")
        

        ##Lancement de la fonction affichageCluster(Dico,fen)
        
        affichageCluster(Tandem,fen2)
        fen2.config(state = Tkinter.DISABLED)
        
        affichageCluster(Segmental,fen3)
        
        ##G�n�ration des fen�tres
        
        


###########################################################################

##Fonction permettant de s�lectionner le chemin d'acc�s au fichier m8################

def ouvrirFichier() :
    #ouverture d'une fenetre de dialogue d'ouverture d'un fichier de type m8 uniquement
    filename = tkFileDialog.askopenfilename(filetypes = [('Fichiers m8', '*.m8')])
      
    if (not filename):
        tkMessageBox.showinfo("Erreur","Fichier non sp�cifi�")    
    else :
    #insertion dans la zone de texte situ�e � c�t� du bouton Parcourir du nom de fichier s�lectionn� apr�s ouverture de la fen�tre de dialogue pr�c�dente
        entrPar.insert(END, filename)
        
##################################################################################

##########Fonction permettant d'acc�der � un lien vers une fiche Uniprot
        
def affichePageWeb():
        
    ##### Copie du code "Creation du dictionnaire SW{idPDB:idSW} (idPDB=id_paramecium ; idSW=lien vers fiche SW)" ######

    ####### Stockage du numero d'accession du gene et sa fiche swiss prot dans un dictionnaire a partir du fichier protAC #######
    ObFichier = open('protAC.txt','r')
    ligne = ObFichier.readline()

    SW={}
    liste1=[]
    liste2=[]
    liste =[]
    liste = ligne.split('\t')   

    idPDB = liste[0]            
    liste1.append(idPDB)        
    idSW = liste[3]             
    idSW=idSW.strip()           
    liste2.append(idSW)         
    SW[idPDB]=idSW                   

    for ligne in ObFichier:     
      liste = ligne.split('\t') 
      if (liste[0]!=''):
        idPDB = liste[0]          
        liste1.append(idPDB)   
        idSW = liste[3]          
        idSW=idSW.strip()        
        liste2.append(idSW)       
        SW[idPDB]=idSW            
    ObFichier.close()
    
########Fin de creation du dictionnaire########

    
#########Recup�ration du lien###############
    try:
        idPDB=str(entrRes.get())#r�cup�re sous forme de cha�ne de carct�res l'identifiant saisi
        idSW=SW[idPDB].replace('_PARTE','')
        webbrowser.open('http://www.uniprot.org/uniprot/'+idSW+'.txt') ##ouvre la page web correspondante
    except:
        tkMessageBox.showinfo("Identifiant de prot�ine non sp�cifi� ou non valide")
#############################################################################################################

###################################Cr�ation de la fen�tre principale#########################################

fen1=Tk()
fen1.config (bg = "#77B5FE")  #d�finition de la couleur de fond

#D�finition d'un premier style de caract�res
fonte = tkFont.Font(fen1)
#family:style de caract�re, slant: normal(ROMAN) ou italique  ,weight: gras ou normal, strikeout: ray� ou non
fonte.config(size=13, family='Andalus', slant=tkFont.ROMAN, weight=tkFont.NORMAL, underline=True, overstrike=False)

#D�finition d'un second style de caract�res
fonte2 = tkFont.Font(fen1)
fonte2.config(size=11, family='Andalus', slant=tkFont.ROMAN, weight=tkFont.NORMAL, underline=False, overstrike=False)

#D�finition d'un troisi�me style de caract�res
fonte3 = tkFont.Font(fen1)
fonte3.config(size=11, family='Consolas', slant=tkFont.ROMAN, weight=tkFont.NORMAL, underline=False, overstrike=False)

fen1.title('Duplications chez Paramecium tetraurelia')


#Cr�ation d'un texte indiquant de s�lectionner un fichier
txt1=Label(fen1,bg = "#77B5FE",text='S�lectionner un fichier de sortie d\'alignement Blastall (format m8)')
txt1.config(font = fonte2)#configuration du style de caract�res fonte2
#zone de texte d'affichage du fichier s�lectionn� par Parcourir
entrPar=Entry(fen1,width=50)
#Cr�ation d'un bouton Parcourir
bou1=Button(fen1,bg='grey',text='Parcourir',width=10, command=ouvrirFichier)
bou1.config(font = fonte3)
#Cr�ation d'un texte "Param�tres d'ajustement"
txt2=Label(fen1,bg = "#77B5FE",text='Param�tres d\'ajustement')
txt2.config(font = fonte)
#Cr�ation d'un texte "E-value"
txt2a=Label(fen1,bg = "#77B5FE",text='E-value')
txt2a.config(font = fonte2)
#Cr�ation d'une zone de texte pour indiquer une e-value
entr2a=Entry(fen1)
#Cr�ation d'un texte "Pourcentage d'identit�"
txt2b=Label(fen1,bg = "#77B5FE",text='Pourcentage d\'identit�')
txt2b.config(font = fonte2)
#Zone de texte pour indiquer un pourcentage d'identit�
entr2b=Entry(fen1)
#Cr�ation d'un bouton "Lancer"
bou2=Button(fen1,text='Lancer',width=10, command=verif_value)
bou2.config(font = fonte3)
#Cr�ation d'un texte titre "Informations"
txtRes=Label(fen1, bg = "#77B5FE", text='Informations d\'une prot�ine')
txtRes.config(font = fonte)
#Cr�ation d'une zone de texte d'insertion d'un id prot�ine
entrRes=Entry(fen1)
#Cr�ation d'un bouton renvoyant � un lien sur une fiche Uniprot
bou3=Button(fen1,text='Fiche Uniprot',width=15, command=affichePageWeb)
bou3.config(font = fonte3)

#Disposition des �l�ments dans la fen�tre
txt1.grid(row = 0,column = 1)
txt2.grid(row = 4,column = 1)
txt2a.grid(row = 5, column = 1, sticky = E) 
txt2b.grid(row = 7, column = 1, sticky = E)
entrPar.grid(row = 2, column =1, sticky = E)
entr2a.grid(row = 5,column = 2, sticky = W)
entr2b.grid(row = 7,column = 2, sticky = W)
bou1.grid(row = 2, column = 2, sticky = W)
bou2.grid(row = 8,column = 2, sticky = W)
txtRes.grid(row = 9, column = 1)
entrRes.grid(row= 10, column = 1, sticky = E)
bou3.grid(row = 10,column = 2, sticky = W)


#Cr�ation d'un menu 'Fichier'
mainmenu = Menu(fen1)  
menuFichier = Menu(mainmenu)
menuFichier.add_command(label="Quitter", command=fen1.quit)
mainmenu.add_cascade(label = "Fichier", menu = menuFichier)
fen1.config(menu = mainmenu)
    
fen1.mainloop()#G�n�ration de la fen�tre principale
