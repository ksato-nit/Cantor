<?php
$input = file_get_contents("extended_mumford_projective.cpp");
$rows = explode("\n", $input);
$output = "";

foreach($rows as $ind => $row){
    $pattern = "/.*mpz_sub\([a-zA-Z0-9].*\, [a-zA-Z0-9].*\, [a-zA-Z0-9].*\).*/";
    $match = preg_match($pattern, $row, $matches);
    if($match != false){
        $start = strpos($row, "mpz_sub");
        $params = substr($row, $start + 8, strpos($row, ")") - $start - 8);
        $params = explode(", ", $params);
        // $op1 = $op2 + $op3.
        if(count($params) == 2){
            var_dump($row);
            var_dump($params);
        }
        $op1 = substr($params[0], 0, strpos($params[0], ".value"));
        $op2 = substr($params[1], 0, strpos($params[1], ".value"));
        $op3 = substr($params[2], 0, strpos($params[2], ".value"));
        $row1 = "mpz_sub(" . $op1 . ".re, " . $op2 . ".re, " . $op3 . ".re);";
        $row2 = "mpz_sub(" . $op1 . ".im, " . $op2 . ".im, " . $op3 . ".im);";
        $output .= $row1 . "\n" . $row2 . "\n";
    }else{
        $output .= $row . "\n";
    }
}

file_put_contents("extended_mumford_projective.cpp", $output);
