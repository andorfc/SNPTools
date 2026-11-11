<?php
/* Runs one SNPTools PHP endpoint as if requested over HTTP (CLI stand-in for php -S).
   Request is passed as JSON in SNPT_REQ: {method, get, post, root, script}. */
$req = json_decode(getenv('SNPT_REQ'), true);
$_SERVER['REQUEST_METHOD'] = $req['method'];
$_GET = $req['get'] ?: array();
$_POST = $req['post'] ?: array();
$_REQUEST = array_merge($_GET, $_POST);
chdir($req['root']);
include $req['script'];
