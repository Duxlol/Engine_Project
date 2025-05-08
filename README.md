## Gequoteerde functionaliteit

V: Werkend  
-: Deels werkend met gekende problemen (onderaan beschreven)  
X: Niet werkend of niet geïmplementeerd  
(blanco): TODO  

|   | Functionaliteit               | Status |
|---|-------------------------------|--------|
| 1 | **2D L-systemen**             | V      |
|   | Met haakjes                   | V      |
|   | Stochastisch                  | X      |
| 2 | Transformaties                |        |
|   | Eye-point                     | X      |
|   | Projectie                     | X      |
| 3 | Platonische Lichamen          | X      |
|   | Kegel en cylinder             | X      |
|   | Bol                           | X      |
|   | Torus                         | X      |
|   | 3D L-systemen                 | X      |
| 4 | Z-buffering (lijnen)          | X      |
| 5 | Triangulatie                  | X      |
|   | Z-buffering (driehoeken)      | X      |
| 6 | 3D fractalen                  | X      |
|   | BuckyBall                     | X      |
|   | Mengerspons                   | X      |
|   | View Frustum                  | X      |
| 7 | Ambient licht                 | X      |
|   | Diffuus licht (oneindig)      | X      |
|   | Diffuus licht (puntbron)      | X      |
|   | Speculair licht               | X      |
| 8 | Schaduw                       | X      |
|   | Texture mapping               | X      |
| 9 | Bollen en cylinders           | X      |
|   | UV-coordinaten                | X      |
|   | Cube mapping                  | X      |

Geïmplementeerde vorm van texture mapping: *Niet geïmplementeerd*

---

## Gekende problemen 
- **Haakjes in L-systemen**: Er kunenn kleine afwijkingen ontstaan bij complexe vertakkingen door fouten bij afrondingen in de hoekberekeningen
- **Schaalafhankelijkheid**: Bij heel grote L-systemen kunnen lijnen buiten het beeld vallen door fixed scaling

---

## Niet-gequoteerde functionaliteit
- Ondersteuning voor aangepaste achtergrondkleuren via INI-configuratie.
- Dynamische schaling van L-systemen om binnen de beeldafmetingen te passen.

---

## Extra functionaliteit, niet in de opgaves beschreven
