# esta gaita tem swag como o caracas, assim que tiver tempo tenho que fazer um tips and tricks de sed e meter uma pagina no cheat como deve de ser!
cat $DATAFILE | (sed -u '/redshift,luminosity_distance,error/'q; sort -n)

# explicacao:
# - imprimir ficheiro para STDOUT
# - pipe para um parentisses com a syntax (command A; command B). esta syntax significa que o STDIN vai ser lido pelo 'command A' ate ele acabar, ficando o 'command B' com o resto
# - aqui mandamos para o sed (-u: unbuffered, ou seja manda logo nao fica a espera de receber tudo) que quando ler a frase dentro das quotes vai sair (q no final da string)
# - o resto, que vao ser apenas os numeros do .csv, vao para o sort (-n para dar sort pelo valor numerico), que vai dar sort ao resto do ficheiro em funcao do 1º numero

# cuidado para nao mandares o output para o mesmo ficheiro pq senao ele vai ficar logo vazio! manda para um ficheiro temporario (no /tmp se quiseres ser roleplayer) e dps da overwrite ao ficheiro original
