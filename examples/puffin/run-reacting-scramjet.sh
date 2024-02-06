#! /bin/bash
puffin-prep --job=reacting-scramjet
puffin --job=reacting-scramjet
puffin-post --job=reacting-scramjet --output=vtk

