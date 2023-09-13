import random
import math
import matplotlib.pyplot as plt
hydrophobic_aa = "FAMILYVWC"

class Amino_acid :

    def __init__(self,numero,type,x=0.0,y=0.0):
        if (type in hydrophobic_aa):
            self.type = "H"
        else:
            self.type = "P"
        self.numero = numero
        self.x = x
        self.y = y

    def __str__(self):
        return f"Il s'agit de l'acide aminé numéro {self.numero} de type {self.type} et de coordonnées [{self.x},{self.y}]"

    def get_type (self):
        return self.type

    def get_numero (self):
        return self.numero

    def get_coord (self):
        return self.x, self.y

    def set_coord (self,new_x,new_y):
        self.x = new_x
        self.y = new_y



class Protein :
    def __init__(self,filename):
       self.sequence_aa = self.lecture_seq(filename)
       self.filename = filename
       self.current_energy = 0
       self.initial_structure()

    def initial_structure(self):
        """ Initialise la structure initiale de la protéine en respectant les mouvements implémentés."""

        l = len(self.sequence_aa)

        # Place le premier acide aminé (n°1) à la position (0, 0)
        self.sequence_aa[0].set_coord(0, 0)

        # Place les acides aminés hydrophobes (H) de manière à ce que les mouvements possibles soient respectés
        for i in range(1, l):
            prev_x, prev_y = self.sequence_aa[i - 1].get_coord()
            if i % 2 == 1:
                # Pour les acides aminés impairs, ils seront placé a coté du précédent
                self.sequence_aa[i].set_coord(prev_x + 1, prev_y)
            else:
                # Pour les acides aminés pairs, ils sont placés au dessus du précédent
                if i % 4 == 0:
                    self.sequence_aa[i].set_coord(prev_x, prev_y + 1)
                else:
                    self.sequence_aa[i].set_coord(prev_x, prev_y - 1)

    def lecture_seq(self,filename):
        """Fonction qui permet de lire un fichier fasta et de traduire les acides aminées
        selon le modèle HP, cette fonction crée également un nouveau fichier au format fasta
        avec les acides aminées modifiées.

        args :
        - filename : nom du fichier fasta """

        prot_formatee = ""
        sequence = []
        with open(filename, "r") as fasta:
            with open(filename[:-6] + '_HP' + filename[-6:], "w") as out:
                for line in fasta:
                    if not line.startswith(">"):
                        for aa in line.strip():
                            if aa in hydrophobic_aa:
                                aa = "H"
                            else:
                                aa = "P"
                            out.write(f'{aa}')
                            prot_formatee = prot_formatee + aa

        for indice, aa in enumerate(prot_formatee):
            sequence.append(Amino_acid(indice + 1, aa, x=indice, y=0))

        return sequence


    def full_seq(self):
        seq = []
        for aa in self.sequence_aa:
            seq.append([aa.get_numero(),aa.get_type(), aa.get_coord()])
        return seq


    def energetic_measure(self):
        """fonction qui calcul l'energie de la protéine en 2D"""
        l = len(self.sequence_aa)

        for i in range(l):
            for j in range(i + 1, l):
                if self.sequence_aa[i].get_type() == 'H' and self.sequence_aa[j].get_type() == 'H':
                    distance = self.calculate_distance(self.sequence_aa[i], self.sequence_aa[j])
                    if self.is_valid_hh_contact(distance):
                        self.current_energy -= 1  # Chaque contact H-H valide contribue -1 à l'énergie totale

        return self.current_energy

    def calculate_distance(self, aa1, aa2):
        """Fonction qui permet de calculer la distance euclidienne entre deux acides aminés """
        x1, y1 = aa1.get_coord()
        x2, y2 = aa2.get_coord()
        distance = math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
        return distance


    def is_valid_hh_contact(self,distance):
        if distance == 1:
            return True
        else:
            return False


    def corner_move(self):
        """Fonction qui permet d'effectuer un corner move si cela est permis"""
        l = len(self.sequence_aa)
        if l < 3:
            return "Attention il n'y a pas assez de résidus dans la séquence"

        # choix d'un résidu au hasard en dehors des extremités
        i = random.randint(1,l-2)

        # verification si les résidus voisins (i -1 et i + 1) partagent un voisin commun
        # voisins en i -1
        neighbor_position_i_moins_1 = set([
            (self.sequence_aa[i - 1].x, self.sequence_aa[i - 1].y + 1),
            (self.sequence_aa[i - 1].x, self.sequence_aa[i - 1].y - 1),
            (self.sequence_aa[i - 1].x + 1, self.sequence_aa[i - 1].y),
            (self.sequence_aa[i - 1].x - 1, self.sequence_aa[i - 1].y)
        ])

        neighbor_position_i_plus_1 = set([
            (self.sequence_aa[i + 1].x, self.sequence_aa[i + 1].y + 1),
            (self.sequence_aa[i + 1].x, self.sequence_aa[i + 1].y - 1),
            (self.sequence_aa[i + 1].x + 1, self.sequence_aa[i + 1].y),
            (self.sequence_aa[i + 1].x - 1, self.sequence_aa[i + 1].y)
        ])

        common_neighbors = neighbor_position_i_moins_1.intersection(neighbor_position_i_plus_1)


        if not common_neighbors:
            return "Attention il n'y a pas assez de résidus dans la séquence"

        # Choix d'un voisin commun
        new_x,new_y = random.choice(list(common_neighbors))

        #Deplacement du résidu i vers sa nouvelle position
        self.sequence_aa[i].set_coord(new_x,new_y)

        #Calcul de l'energie post-movement
        self.current_energy = self.energetic_measure()

        return "Mouvement possible"


    def end_move(self):
        """Fonction qui permet d'effectuer une end move si cela est permis'"""

        n = len(self.sequence_aa)

        if (n < 2):
            return ("Attention il n'y a pas assez de résidus dans la séquence")

        # Choix du premier ou dernier résidu (1 ou n)
        end_to_move = random.choice([1,n])

        if (end_to_move == 1):
            # Cas ou le premier résidu est choisi
            new_x,new_y = self.sequence_aa[0].x, self.sequence_aa[0].y

            while (new_x,new_y) == (self.sequence_aa[0].x, self.sequence_aa[0].y):
                # Choix d'une position adjacente au premier residus'

                neighbors = [(self.sequence_aa[1].x, self.sequence_aa[1].y + 1),
                             (self.sequence_aa[1].x, self.sequence_aa[1].y - 1),
                             (self.sequence_aa[1].x + 1, self.sequence_aa[1].y),
                             (self.sequence_aa[1].x - 1, self.sequence_aa[1].y)]
                # Choix aléatoire d'une position adjacente au premier residus
                new_x, new_y = random.choice(neighbors)

        else:
            # Cas ou le dernier résidu est choisi
            new_x, new_y = self.sequence_aa[-1].x, self.sequence_aa[-1].y

            while (new_x,new_y) == (self.sequence_aa[-1].x, self.sequence_aa[-1].y):
                # Choix d'une position adjacente au dernier residus'

                neighbors = [(self.sequence_aa[-2].x, self.sequence_aa[-2].y + 1),
                             (self.sequence_aa[-2].x, self.sequence_aa[-2].y - 1),
                             (self.sequence_aa[-2].x + 1, self.sequence_aa[-2].y),
                             (self.sequence_aa[-2].x - 1, self.sequence_aa[-2].y)]
                # Choix aléatoire d'une position adjacente au dernier residus
                new_x, new_y = random.choice(neighbors)

            self.sequence_aa[-1].set_coord(new_x, new_y)

        #Calcul de l'energie post-movement
        self.current_energy = self.energetic_measure()

        return "Mouvement possible"


    def crankshaft_move(self):
        """Fonction qui permet d'effectuer un crankshaft move si cela est permis'"""

        n = len(self.sequence_aa)
        if (n < 4):
            return "La séquence est trop courte pour pouvoir effectuer un crankshaft move"

        # Choix d'un résidu au hasard en dehors des extremités
        i = random.randint(1,n-3)

        # Vérfie si les positions i' et i' + 1 sont vides
        position_i = (self.sequence_aa[i].x, self.sequence_aa[i].y)
        position_i_prime = (self.sequence_aa[i+2].x, self.sequence_aa[i + 2].y)


        if position_i_prime not in [(position_i[0] + 1, position_i[1]),
                                    (position_i[0] - 1, position_i[1]),
                                    (position_i[0], position_i[1] + 1),
                                    (position_i[0], position_i[1] - 1)]:
            # Les position i' ne sont pas des voisins de i
            return "Mouvement crankshaft impossible"

        # Rotation de 180 degrés de la structure en U
        self.sequence_aa[i].set_coord(self.sequence_aa[i + 2].x, self.sequence_aa[i + 2].y)
        self.sequence_aa[i + 1].coord(self.sequence_aa[i + 1].x, self.sequence_aa[i + 1].y)
        self.sequence_aa[i + 2].set_coord(position_i[0], position_i[1])


        # Calcul de l'energie post-movement
        self.current_energy = self.energetic_measure()

        return "Mouvement possible"

    def minimize_energy(self, max_iterations=100):
        """Minimise l'énergie globale de la protéine en effectuant des mouvements aléatoires par acide aminé."""
        best_energy = self.energetic_measure()  # Énergie initiale
        best_sequence = list(self.sequence_aa)  # Configuration initiale

        for i in range(max_iterations):
            # Sauvegarde de l'état actuel de la séquence
            current_sequence = list(self.sequence_aa)
            current_energy = self.current_energy

            # Parcourir chaque acide aminé et appliquer un mouvement aléatoire
            for aa in self.sequence_aa:
                # Choix aléatoire d'un mouvement parmi ceux disponibles pour cet acide aminé
                available_moves = [self.corner_move, self.crankshaft_move, self.end_move]
                move_choice = random.choice(available_moves)

                # Appliquer le mouvement choisi à l'acide aminé
                result = move_choice()

                if result != "Mouvement possible":
                    # Si le mouvement est impossible, restaurer l'état précédent
                    self.sequence_aa = list(current_sequence)
                    self.current_energy = current_energy

            # Calcul de la nouvelle énergie après chaque itération
            new_energy = self.energetic_measure()

            # Mettre à jour la séquence si l'énergie est meilleure ou avec une certaine probabilité si l'énergie est pire
            if new_energy < best_energy or random.random() < math.exp(-(new_energy - best_energy)):
                best_energy = new_energy
                best_sequence = list(self.sequence_aa)

        # Mettre à jour la séquence avec la meilleure configuration trouvée
        self.sequence_aa = list(best_sequence)
        self.current_energy = best_energy

    def visualize_sequence_2d(self):
        plt.figure()
        x = [aa.get_coord()[0] for aa in self.sequence_aa]
        y = [aa.get_coord()[1] for aa in self.sequence_aa]


        plt.plot(x, y, 'o', markersize=10, markerfacecolor='lightgray', markeredgecolor='black')

        # Ajout de lignes entre les acides aminés
        for i in range(len(self.sequence_aa) - 1):
            x_coords = [x[i], x[i + 1]]
            y_coords = [y[i], y[i + 1]]
            plt.plot(x_coords, y_coords, '-', color='black')

        plt.title("Représentation 2D de la conformation optimisée de la protéine P02776")
        plt.xlabel('Position X')
        plt.ylabel('Position Y')
        plt.grid(True)

        # Affichage de la grille
        plt.show()

if __name__ == "__main__":
    protein = Protein("P02776.fasta")
    print(protein.full_seq())
    print(protein.energetic_measure())
    protein.initial_structure()
    protein.visualize_sequence_2d()
    protein.minimize_energy()
    protein.visualize_sequence_2d()
    print(f" l'energie minimale est {protein.current_energy}")







