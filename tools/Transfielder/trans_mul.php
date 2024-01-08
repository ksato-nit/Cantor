<?php
$input = file_get_contents("extended_mumford_projective.cpp");
$rows = explode("\n", $input);
$output = "";

foreach($rows as $ind => $row){
    $pattern = "/.*mpz_mul\([a-zA-Z0-9].*\, [a-zA-Z0-9].*\, [a-zA-Z0-9].*\).*/";
    $match = preg_match($pattern, $row, $matches);
    if($match != false){
        $start = strpos($row, "mpz_mul");
        $params = substr($row, $start + 8, strpos($row, ")") - $start - 8);
        $params = explode(", ", $params);
        // $op1 = $op2 + $op3.
        $op1 = substr($params[0], 0, strpos($params[0], ".value"));
        $op2 = substr($params[1], 0, strpos($params[1], ".value"));
        $op3 = substr($params[2], 0, strpos($params[2], ".value"));
        if($op1 == $op2 || $op1 == $op3){
            // このとき中間変数を用意する必要がある
            $row1 = "mpz_mul(tempd.re, " . $op2 . ".re, " . $op3 . ".re);";
            $row2 = "mpz_mod(tempd.re, tempd.re, ExtendedNumber::CHARA);" ;
            $row3 = "mpz_mul(temp, " . $op2 . ".im, " . $op3 . ".im);";
            $row4 = "mpz_mod(temp, temp, ExtendedNumber::CHARA);" ;
            $row5 = "mpz_sub(tempd.re, tempd.re, temp);";

            $row6 = "mpz_mul(tempd.im, " . $op2 . ".re, " . $op3 . ".im);";
            $row7 = "mpz_mod(tempd.im, tempd.im, ExtendedNumber::CHARA);" ;
            $row8 = "mpz_mul(temp, " . $op2 . ".im, " . $op3 . ".re);";
            $row9 = "mpz_mod(temp, temp, ExtendedNumber::CHARA);" ;
            $row10 = "mpz_add(tempd.im, tempd.im, temp);";

            $row11 = "mpz_set(" . $op1 . ".re, tempd.re);";
            $row12 = "mpz_set(" . $op1 . ".im, tempd.im);";
        }else{
            $row1 = "mpz_mul(" . $op1 . ".re, " . $op2 . ".re, " . $op3 . ".re);";
            $row2 = "mpz_mod(" . $op1 . ".re, " . $op1 . ".re, ExtendedNumber::CHARA);";
            $row3 = "mpz_mul(temp, " . $op2 . ".im, " . $op3 . ".im);";
            $row4 = "mpz_mod(temp, temp, ExtendedNumber::CHARA);" ;
            $row5 = "mpz_sub(" . $op1 . ".re, " . $op1 . ".re, temp);";

            $row6 = "mpz_mul(" . $op1 . ".im, " . $op2 . ".re, " . $op3 . ".im);";
            $row7 = "mpz_mod(" . $op1 . ".im, " . $op1 . ".im, ExtendedNumber::CHARA);";
            $row8 = "mpz_mul(temp, " . $op2 . ".im, " . $op3 . ".re);";
            $row9 = "mpz_mod(temp, temp, ExtendedNumber::CHARA);" ;
            $row10 = "mpz_add(" . $op1 . ".im, " . $op1 . ".im, temp);";

            $row11 = "";
            $row12 = "";
        }


        $output .= $row1 . "\n" . $row2 . "\n" . $row3 . "\n" . $row4 . "\n" . $row5 . "\n" . $row6 . "\n" . $row7 . "\n" . $row8 . "\n" . $row9 . "\n" . $row10 . "\n" . $row11 . "\n" . $row12 . "\n";
    }else{
        $output .= $row . "\n";
    }
}

file_put_contents("extended_mumford_projective.cpp", $output);
