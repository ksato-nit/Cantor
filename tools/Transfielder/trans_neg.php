<?php
$input = file_get_contents("extended_mumford_projective.cpp");
$rows = explode("\n", $input);
$output = "";

foreach($rows as $ind => $row){
    $pattern = "/.*mpz_neg\([a-zA-Z0-9].*\, [a-zA-Z0-9].*\).*/";
    $match = preg_match($pattern, $row, $matches);
    if($match != false){
        $start = strpos($row, "mpz_neg");
        $params = substr($row, $start + 8, strpos($row, ")") - $start - 8);
        $params = explode(", ", $params);
        $op1 = substr($params[0], 0, strpos($params[0], ".value"));
        $op2 = substr($params[1], 0, strpos($params[1], ".value"));
        $row1 = "mpz_neg(" . $op1 . ".re, " . $op2 . ".re);";
        $row2 = "mpz_neg(" . $op1 . ".im, " . $op2 . ".im);";
        $output .= $row1 . "\n" . $row2 . "\n";
    }else{
        $output .= $row . "\n";
    }
}

file_put_contents("extended_mumford_projective.cpp", $output);
